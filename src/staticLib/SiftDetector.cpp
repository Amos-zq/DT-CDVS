#include "Defines.h"
#include "SiftDetector.h"

//math operation
float vl_mod_2pi_f (float x) {
    while (x > (float)(CDVS_PI2)) x -= (float) (CDVS_PI2) ;
    while (x < 0.0F) x += (float) (CDVS_PI2);
    return x ;
}
float vl_fast_atan2_f (float y, float x) {
    float angle, r ;
    float abs_y    = fabsf (y) + VL_EPSILON_F ;

    if (x >= 0) {
        r = (x - abs_y) / (x + abs_y) ;
        angle = VL_QUTER_PI;
    } else {
        r = (x + abs_y) / (abs_y - x) ;
        angle = VL_THR_QUTER_PI;
    }
    angle += (0.1821F*r*r - 0.9675F) * r ;
    return (y < 0) ? CDVS_PI2 - angle : angle ;
}
long int vl_floor_f (float x) {
    long int xi = (long int) x ;
    if (x >= 0 || (float) xi == x) return xi ;
    else return xi - 1 ;
}
float vl_fast_resqrt_f (float x) {
    float xhalf = 0.5f * x;
    int i = *(int*)&x;              // get bits for floating VALUE
    i = 0x5f375a86- (i >> 1);       // gives initial guess y0
    x = *(float*)&i;                // convert bits BACK to float
    x = x * (1.5f - xhalf * x * x); // Newton step, repeating increases accuracy
    return x;
}
float vl_fast_sqrt_f (float x) {
    return (x < 1e-8) ? 0 : x * vl_fast_resqrt_f (x) ;
}
long int vl_floor_d (float x) {
    long int xi = (long int) x ;
    if (x >= 0 || (float) xi == x) return xi ;
    else return xi - 1 ;
}
float fast_expn (float x) {
    int i ;
    if (x > EXPN_MAX) return 0.0 ;
    x *= 10.24f;    //EXPN_SZ / EXPN_MAX ;
    i = (int)x;
    return expn_tab[i] + (x - i) * (expn_tab[i + 1] - expn_tab[i]) ;
}

//for both BFlogDetector and FlogDetector
FilterUnit* init_block() {
    FilterUnit * p = new FilterUnit;
    if(!p)  return NULL;
    int fblock_width    = BLOCK_WIDTH + MAX_FILTER_WIDTH - 1;       //the block filtering size
    int imageSize       = fblock_width * fblock_width;				//size of spacial domain
    int fblock_size     = fblock_width * (fblock_width / 2 + 1);	//size of frequency domain

    p->inmat        = (float*) fftwf_malloc (sizeof(float) * imageSize) ;
    p->fblock       = (fftwf_complex*) fftwf_malloc (sizeof(fftwf_complex) * fblock_size) ;
    p->fblock_flog  = (fftwf_complex*) fftwf_malloc (sizeof(fftwf_complex) * fblock_size) ;

    for(int i = 0; i < 2; i++) {
        p->fftplan[i]  = fftwf_plan_dft_r2c_2d(fblock_width, fblock_width, p->inmat,        p->fblock,  FFTW_ESTIMATE);
        p->ifftplan[i] = fftwf_plan_dft_c2r_2d(fblock_width, fblock_width, p->fblock_flog,  p->inmat,   FFTW_ESTIMATE);
    }
    return p;
}
void release_block(FilterUnit * p) {
    if(!p) return;
    if (p->inmat)       fftwf_free (p->inmat);
    if (p->fblock)      fftwf_free(p->fblock);
    if (p->fblock_flog) fftwf_free(p->fblock_flog);

    for (int i = 0; i < 2; ++i) {
        if (p->fftplan[i])  fftwf_destroy_plan(p->fftplan[i]);
        if (p->ifftplan[i]) fftwf_destroy_plan(p->ifftplan[i]);
    }
    delete p;
}

void assign_memory (SiftFilt * f, int width, int height) {
    for (int s = 0; s < f->S + 2; ++s) {
        float sigma = 1.6 * pow (2.0, (s) * 1.0 / f->S);
        f->normFactors[s] = 1.0 * sigma * sigma / 160000.0;		// adjust for scale
    }

    //allocate memory for Gaussian block buffer
    f -> block_row_n = (width / BLOCK_WIDTH + (width % BLOCK_WIDTH == 0 ? 0 : 1));
    f -> block_col_n = (height / BLOCK_WIDTH + (height % BLOCK_WIDTH == 0 ? 0 : 1));
    f -> block_row_remain_width = width + BLOCK_WIDTH - f->block_row_n * BLOCK_WIDTH;
    f -> block_col_remain_height = height + BLOCK_WIDTH - f->block_col_n * BLOCK_WIDTH;
    f->block_buffer_len = f->block_row_n;
    if(f->block_col_n > 1) {
        f->block_buffer_len++;
        if(f->block_row_n > 1) {
            f->block_buffer_len++;
        }
    }
    if(f->block_buffer) free(f->block_buffer);
    f -> block_buffer	= (gs_type *) malloc (sizeof(gs_type) * f->S * BLOCK_SIZE * (f->block_buffer_len));
    for(int s = 0; s < f->block_row_n * f->block_col_n; s++) f->block_map[s] = (s % f->block_buffer_len) * BLOCK_SIZE;
}

node get_filter(SiftFilt const *f, node filter, int s) {
    node res;
    res.pos = filter.pos + s - f->s_min;
    res.mask_pos = filter.mask_pos + s - f->s_min;
    res.data = filter.data + res.pos[0];
    res.mask = filter.mask + res.mask_pos[0];
    return res;
}
node get_filters(int o_cur, int type) {
    node res;
    switch (o_cur) {
    case 0:
        res.data = type ? g_fgaussian_data0:g_flog_data0;
        res.mask = type ? g_fgaussian_mask0:g_flog_mask0;
        res.mask_pos = type ? g_fgaussian_maskpos0:g_flog_maskpos0;
        res.pos = type ? g_fgaussian_pos0:g_flog_pos0;
        break;
    case 1:
        res.data = type ? g_fgaussian_data1:g_flog_data1;
        res.mask = type ? g_fgaussian_mask1:g_flog_mask1;
        res.mask_pos = type ? g_fgaussian_maskpos1:g_flog_maskpos1;
        res.pos = type ? g_fgaussian_pos1:g_flog_pos1;
        break;
    default:
        break;
    }
    return res;
}

void complex_mutil(fftwf_complex *srcA,const node srcB,fftwf_complex *dest,size_t size, int fwidth, float factor) {
    int j = 0, k = 0, h = 0, ii, jj;
    int locate_1, locate_2, locate_3, locate_4;
    int pos1 = 0, pos2 = (fwidth + fwidth - 2) * fwidth;
    memset(dest, 0, sizeof(fftwf_complex)*size);
    for(ii = 0; ii < fwidth; ii++) {
        locate_1 = locate_2 = pos1 + ii;
        locate_3 = locate_4 = pos2 + ii;
        pos1 += fwidth;
        pos2 -= fwidth;
        for(jj = ii; jj < fwidth; jj++) {
            if(srcB.mask[k] & bits[h]) {
                float v = (float)srcB.data[j] * factor;
                dest[locate_1][0] = srcA[locate_1][0] * v;
                dest[locate_1][1] = srcA[locate_1][1] * v;

                dest[locate_2][0] = srcA[locate_2][0] * v;
                dest[locate_2][1] = srcA[locate_2][1] * v;

                if(ii > 0 && ii < fwidth - 1) {
                    dest[locate_3][0] = srcA[locate_3][0] * v;
                    dest[locate_3][1] = srcA[locate_3][1] * v;
                }
                if(jj > 0 && jj < fwidth - 1) {
                    dest[locate_4][0] = srcA[locate_4][0] * v;
                    dest[locate_4][1] = srcA[locate_4][1] * v;
                }
                j++;
            }
            h++;
            if(h == 8) {
                h = 0;
                k++;
            }
            locate_1++;
            locate_2 += fwidth;
            locate_3++;
            locate_4 -= fwidth;
        }
    }
}
float normalize_histogram(float * p) {
    float  norm = 0.0 ;
    for (int i = 0; i < 128; i++)
        norm += p[i] * p[i] ;
    norm = sqrt (norm) + VL_EPSILON_F ;

    if(norm > 0) {
        for (int i = 0; i < 128; i++) {
            p[i] /= norm;
        }
    }
    return norm;
}

//for BFlogDetector
SiftFilt * sift_new (int width, int height, int noctaves, int nlevels, int o_min, unsigned char *data) {
    if(!data)   return NULL;
    SiftFilt *f = (SiftFilt*)malloc (sizeof(SiftFilt)) ;
    if(!f)      return NULL;
    int  LoG_size, Gradient_Size, Input_Size, s;

    if (noctaves < 0) {
        noctaves = (int)VL_MAX (floor (log2 ((float)VL_MIN(width, height))) - o_min - 3, 1);
    }
    f-> width   = width ;
    f-> height  = height ;
    f-> O       = noctaves ;
    f-> S       = nlevels ;
    f-> o_min   = o_min ;
    f-> s_min   = -1 ;
    f-> s_max   = nlevels + 1;
    f-> o_cur   = o_min ;
    f-> sigmak  = pow (2.0, 1.0 / nlevels) ;
    f-> sigma0  = 1.6 * f->sigmak ;

    f-> octave_width  = 0 ;
    f-> octave_height = 0 ;
    f-> block_width = BLOCK_WIDTH;
    f-> maxfilter_width = MAX_FILTER_WIDTH;

    Input_Size = width * height * sizeof(input_type) / 4;
    LoG_size = (BLOCK_WIDTH + 4) * (BLOCK_WIDTH + 4) * (nlevels + 2) * sizeof(element_type);
    Gradient_Size = 120 * 120 * 2 * sizeof(sift_pix);
    f-> data0           = data;
    f-> data            = (input_type*) malloc(Input_Size);
    f-> buffer			= (char *) malloc (VL_MAX(LoG_size, Gradient_Size));
    f-> block_buffer    = NULL;
    f-> log_response    = (element_type *)f->buffer;
    f-> grad			= (sift_pix *)f->buffer;
    f-> grad_map         = NULL;

    #pragma omp critical
    {
        f-> filter	= init_block();
    }
    memset(&(f->fgaussian_filter), 0, sizeof(node));
    memset(&(f->flog_filter), 0, sizeof(node));

    f-> peak_thresh = 0.5;
    f-> edge_thresh = 10.0;
    f-> norm_thresh = 0.0 ;
    f-> magnif      = 3.0 ;
    f-> windowSize  = NBP / 2 ;

    for (s = 0; s < nlevels + 2; ++s) {
        float sigma = 1.6 * pow (2.0, (s) * 1.0 / nlevels);
        f->normFactors[s] = 1.0 * sigma * sigma / 160000.0;		// adjust for scale
    }

    f-> keys     = NULL ;
    f-> nkeys    = 0 ;
    f-> keys_res = 0 ;
    return f ;
}
void sift_delete(SiftFilt* f) {
    if (f) {
        if (f->data) free(f->data);
        if (f->buffer) free(f->buffer);
        if (f->block_buffer) free(f->block_buffer);
#pragma omp critical
        {
            if(f->filter) release_block(f->filter);

        }
        if (f->keys) free (f->keys) ;
        free (f) ;
    }
}
int keypoint_detect (SiftFilt *f, FeatureList & featurelist, float & filter_time, float & extrema_time, float &descript_time) {
    int x, y, s, i, j, ii, jj, iii, jjj;
    int out_cnt = 0;
    /* shortcuts */
    int width           = f-> width ;
    int height          = f-> height;
    int o_min           = f-> o_min ;
    int s_min           = f-> s_min ;
    int s_max           = f-> s_max ;
    int o_cur           = f-> o_cur  ;
    int S				= f-> S;
    int O				= f-> O;
    int w = f-> octave_width  = VL_SHIFT_LEFT(f->width,  - o_cur) ;	//width of current octave
    int h = f-> octave_height = VL_SHIFT_LEFT(f->height, - o_cur) ;	//height of current octave
    int w2 = w / 2;
    int h2 = h / 2;

    float xper  = pow (2.0, o_cur) ;
    float * pdata		= NULL;
    float * pdata2		= NULL;
    gs_type * p		    = NULL;
    gs_type * pt		= NULL;
    input_type * p2		= NULL;
    input_type * pt2	= NULL;
    sift_pix * pgrad	= NULL;
    sift_pix * pgrad2	= NULL;
    element_type *pLoG	= NULL;
    element_type *p2LoG = NULL;
    input_type *data	= f->data;
    unsigned char *data0= f->data0;
    unsigned char *p0	= NULL;
    element_type *log	= f->log_response;

    int h_mode, v_mode;
    int gradient_w, gradient_h, gradient_sx, gradient_sy, gradient_size;
    int extrema_w, extrema_h, extrema_sx, extrema_sy, extrema_size;
    int R_Gradient = MAX_FILTER_WIDTH / 2;
    int R_Extrema = 2;
    int st_key;
    SiftKeypoint *k ;
    int out_total = 0, block_cnt = 0;
    int block_width, block_height;

    assign_memory(f, w, h);
    s_max           = f-> s_max ;
    S				= f-> S;
    int ocur        = o_cur > 0 ? 1 : o_cur;
    int blockw      = f->block_width;		//block width without padding
    int blocksize   = blockw * blockw;
    int padWidth    = f->maxfilter_width / 2;		//padding width
    int fblockw     = blockw + padWidth + padWidth;
    int fblocksize  = fblockw * (fblockw / 2 + 1);

    float fft_scale_log_factor  = 1.0 / (fblockw * fblockw) / FACTOR_LOG;
    float fft_scale_factor      = 1.0 / (fblockw * fblockw) / FACTOR_;


    float *inmat                = f->filter->inmat;
    fftwf_complex *fblock       = f->filter->fblock;
    fftwf_complex *fblock_flog  = f->filter->fblock_flog;
    fftwf_plan  forward_plan    = f->filter->fftplan[ocur];
    fftwf_plan  backward_plan   = f->filter->ifftplan[ocur];
    f->flog_filter      =  get_filters(ocur, 0);
    f->fgaussian_filter =  get_filters(ocur, 1);

    for(y = 0; y < h; y += blockw) {
        v_mode = 0;									    //locate the boundary of the block
        gradient_sy = padWidth - R_Gradient;			//the Y-offset of the gradient patch to the padded block
        gradient_h = blockw + R_Gradient + R_Gradient;  //the height of the gradient patch, with R padding, R <= padWidth
        extrema_sy = padWidth - R_Extrema;
        extrema_h = blockw + R_Extrema + R_Extrema;	    //the height of the extrema patch
        block_height = blockw;
        if(y == 0)	 {							        //the condition of up bound
            v_mode |= 1;
            gradient_sy += R_Gradient;		            //no padding
            gradient_h -= R_Gradient;
            extrema_sy += R_Extrema;
            extrema_h -= R_Extrema;
        }
        if(y + blockw >= h) {
            v_mode |= 2;
            gradient_h = gradient_h - blockw - R_Gradient + h - y;
            extrema_h = extrema_h - blockw - R_Extrema + h - y;
            block_height = h - y;
        }
        for(x = 0; x < w; x += blockw)  {
            h_mode = 0;
            gradient_sx = padWidth - R_Gradient;
            gradient_w = blockw + R_Gradient + R_Gradient;
            extrema_sx = padWidth - R_Extrema;
            extrema_w = blockw + R_Extrema + R_Extrema;
            block_width = blockw;
            if(x == 0) {
                h_mode |= 1;
                gradient_sx += R_Gradient;
                gradient_w -= R_Gradient;
                extrema_sx += R_Extrema;
                extrema_w -= R_Extrema;
            }
            if(x + blockw >= w) {
                h_mode |= 2;
                gradient_w = gradient_w - blockw - R_Gradient + w - x;
                extrema_w = extrema_w - blockw - R_Extrema + w - x;
                block_width = w - x;
            }
            extrema_size = extrema_w * extrema_h;
            gradient_size = gradient_w * gradient_h;

            memset(inmat, 0, sizeof(float) * fblockw * fblockw);
            for(j = y - padWidth, jjj = 0; j < y + blockw + padWidth; j++, jjj++) {
                jj = j;
                //mirror padding of the boundary
                if(jj < 0) {jj = -jj-1;if(jj >= h) continue;}
                if(jj >= h) {jj = h - jj + h - 1;if(jj < 0) continue;}

                if(o_cur == 0) {
                    p0 = data0 + jj * w;
                    pdata = inmat + jjj * fblockw;
                    for(i = x - padWidth, iii = 0; i < x + blockw + padWidth; i++, iii++) {
                        ii = i;
                        //mirror padding of the boundary
                        if(ii < 0) {ii = -ii-1;if(ii >= w) continue;}
                        if(ii >= w){ii =w - ii + w - 1;if(ii < 0) continue;}
                        *(pdata+iii) = (float)*(p0 + ii);
                    }
                } else {
                    p2 = data + jj * w;
                    pdata = inmat + jjj * fblockw;
                    for(i = x - padWidth, iii = 0; i < x + blockw + padWidth; i++, iii++) {
                        ii = i;
                        //mirror padding of the boundary
                        if(ii < 0) {ii = -ii-1;if(ii >= w) continue;}
                        if(ii >= w){ii =w - ii + w - 1;if(ii < 0) continue;}
                        *(pdata+iii) = (float)*(p2 + ii) * FACTOR_;
                    }
                }

            }
            fftwf_execute(forward_plan);
            /* Filtering the block using frequency log filter */
            for(s = s_min; s < s_max ; ++s) {
                complex_mutil(fblock, get_filter(f, f->flog_filter, s), fblock_flog, fblocksize, fblockw / 2 + 1, f->normFactors[s-s_min]);
                fftwf_execute(backward_plan);

                pLoG = f->log_response + (s - s_min) * extrema_size;
                pdata = inmat + extrema_sy * fblockw;
                for(j = 0, jjj = extrema_sy; j < extrema_h; j++, jjj++) {
                    for(i = 0, iii = extrema_sx; i < extrema_w; i++, iii++) {
                        float t = 0;
                        if(iii - 1 >= 0) t += *(pdata+iii) - *(pdata+iii-1);
                        if(iii + 1 < fblockw) t += *(pdata+iii) - *(pdata+iii+1);
                        if(jjj - 1 >= 0) t += *(pdata+iii) - *(pdata+iii-fblockw);
                        if(jjj + 1 < fblockw) t += *(pdata+iii) - *(pdata+iii+fblockw);
                        *(pLoG+i) = (element_type)(t * fft_scale_log_factor);
                    }
                    pLoG += extrema_w;
                    pdata += fblockw;
                }
            }

            //detect block
            st_key = keypoint_extrema_block(f, extrema_w, extrema_h, v_mode, h_mode, R_Extrema);

            st_key = f->nkeys - st_key;     //the keypoint detected of this octave within this block
            /* Filtering the block using gaussian filters */
            for (s = s_min; s < s_min+S; ++s) {
                complex_mutil(fblock, get_filter(f, f->fgaussian_filter, s), fblock_flog, fblocksize, fblockw / 2 + 1, FACTOR);
                fftwf_execute(backward_plan);

                //set block and save the data to the block_buffer;
                p = f->block_buffer + (s - s_min) * f->block_buffer_len * BLOCK_SIZE + f->block_map[block_cnt];
                pdata = inmat + padWidth * fblockw;
                for(j = 0, jjj = padWidth; j < block_height; j++, jjj++) {
                    for(i = 0, iii = padWidth; i < block_width; i++, iii++)
                        *(p+i) = ((gs_type)((int)((*(pdata+iii)) * fft_scale_factor + 0.5)));
                    p += blockw;
                    pdata += fblockw;
                }

                //Copy the last scale to the data and down-sampling
                if(o_cur < f->O - 1 && s == s_min + S - 1) {
                    pdata = inmat + padWidth * fblockw;
                    int W = w;
                    if(o_cur == 0) W = w2;
                    pt2 = data + y / 2 * W;
                    for(j = y / 2, jjj = padWidth; jjj < blockw + padWidth && j < h2; j++, jjj += 2) {
                        for(i = x / 2, iii = padWidth; iii < blockw + padWidth && i < w2; i++, iii += 2)
                            *(pt2+i) = ((input_type)((int)((*(pdata+iii)) * fft_scale_factor + 0.5)));
                        pdata += fblockw + fblockw;
                        pt2 += W;
                    }
                }

                for(i = st_key; i < f->nkeys; i++) {
                    int w;
                    if(f->keys[i].is != s+1)    continue;
                    k = f->keys + i;
                    k->block_id = block_cnt;
                    k->ix += extrema_sx - padWidth;
                    k->iy += extrema_sy - padWidth;
                    k->x += extrema_sx - padWidth;
                    k->y += extrema_sy - padWidth;
                    k->mode = 0;

                    w = (int)(k->sigma * 10.6066 + 0.5);

                    if(!(h_mode & 1) && k->ix - w < 0)
                        k->mode |= MODE_LEFT;
                    if(!(h_mode & 2) && k->ix + w >= block_width)
                        k->mode |= MODE_RIGHT;
                    if(!(v_mode & 1) && k->iy - w < 0)
                        k->mode |= MODE_TOP;
                    if(!(v_mode & 2) && k->iy + w >= block_height)
                        k->mode |= MODE_BOTTOM;
                }
            }
            keypoint_block_desc(f, block_cnt, block_width, block_height, featurelist, xper);//*/
            block_cnt++;
        }//end of block x
    }//end of block y
    if(o_cur > 0) {
        p2 = data;
        pt2 = data;
        for(j = 0; j < h2; j++) {
            for(i = 0; i < w2; i++)
                *(pt2+i) = *(p2+i);
            p2 += w;
            pt2 += w2;
        }
    }
    return 0 ;
}
int keypoint_extrema_block (SiftFilt * f, int w, int h, int v_mode, int h_mode, int pad) {
    element_type* log   = f-> log_response ;
    element_type * pt = NULL;
    int          s_min = f-> s_min ;
    int          s_max = f-> s_max ;
    float       te    = f-> edge_thresh ;
    float       tp    = f-> peak_thresh / FACTOR_LOG;

    int offset;

    int const    xo    = 1 ;      /* x-stride */
    int const    yo    = w ;      /* y-stride */
    int const    so    = w * h ;  /* s-stride */

    int x, y, s, i, ii, jj ;
    int l, r, u, d;
    sift_pix v ;
    SiftKeypoint *k ;
    int st_key;

    /* clear current list */
    //f-> nkeys = 0 ;
    st_key = f->nkeys;

    if(f->o_cur == 0) offset = 6;
    else if(f->o_cur == 1) offset = 3;
    else offset = 2;

    if(h_mode & 1) l = offset; else l = pad;
    if(h_mode & 2) r = offset; else r = pad;
    if(v_mode & 1) u = offset; else u = pad;
    if(v_mode & 2) d = offset; else d = pad;

    /* -----------------------------------------------------------------
    *                                          Find local maxima of DoG
    * -------------------------------------------------------------- */

    /* start from dog [1,1,s_min+1] */

    pt = log + l * xo + u * yo + so ;

    for(s = s_min + 1 ; s <= s_max - 2 ; ++s)
    {
        for(y = u ; y < h - d; ++y) {
            for(x = l ; x < w - r ; ++x) {
                v = *pt ;

#define CHECK_NEIGHBORS(CMP,SGN)                    \
    ( v CMP ## = SGN 0.8 * tp &&                \
    v CMP *(pt + xo) &&                       \
    v CMP *(pt - xo) &&                       \
    v CMP *(pt + so) &&                       \
    v CMP *(pt - so) &&                       \
    v CMP *(pt + yo) &&                       \
    v CMP *(pt - yo) &&                       \
    \
    v CMP *(pt + yo + xo) &&                  \
    v CMP *(pt + yo - xo) &&                  \
    v CMP *(pt - yo + xo) &&                  \
    v CMP *(pt - yo - xo) &&                  \
    \
    v CMP *(pt + xo      + so) &&             \
    v CMP *(pt - xo      + so) &&             \
    v CMP *(pt + yo      + so) &&             \
    v CMP *(pt - yo      + so) &&             \
    v CMP *(pt + xo      - so) &&             \
    v CMP *(pt - xo      - so) &&             \
    v CMP *(pt + yo      - so) &&             \
    v CMP *(pt - yo      - so) &&             \
    v CMP *(pt + yo + xo + so) &&             \
    v CMP *(pt + yo - xo + so) &&             \
    v CMP *(pt - yo + xo + so) &&             \
    v CMP *(pt - yo - xo + so) &&             \
    v CMP *(pt + yo + xo - so) &&             \
    v CMP *(pt + yo - xo - so) &&             \
    v CMP *(pt - yo + xo - so) &&             \
    v CMP *(pt - yo - xo - so) )//*/

                if (CHECK_NEIGHBORS(>,+) ||
                    CHECK_NEIGHBORS(<,-) ) {

                        /* make room for more keypoints */
                        if (f->nkeys >= f->keys_res) {
                            f->keys_res += 200 ;
                            if (f->keys) {
                                f->keys = (SiftKeypoint *)realloc (f->keys,
                                    f->keys_res *
                                    sizeof(SiftKeypoint)) ;
                            } else {
                                f->keys = (SiftKeypoint *)malloc (f->keys_res *
                                    sizeof(SiftKeypoint)) ;
                            }
                        }

                        k = f->keys + (f->nkeys ++) ;

                        k-> ix = x ;
                        k-> iy = y ;
                        k-> is = s ;

                        k-> peak = (float)fabsf(v);

                }
                pt += 1 ;
            }
            pt += l + r ;
        }
        pt += (u + d) * yo ;
    }

    /* -----------------------------------------------------------------
    *                                               Refine local maxima
    * -------------------------------------------------------------- */

    /* this pointer is used to write the keypoints back */
    k = f->keys + st_key;

    for (i = st_key ; i < f->nkeys ; ++i) {

        int x = f-> keys [i] .ix ;
        int y = f-> keys [i] .iy ;
        int s = f-> keys [i]. is ;

        float Dx=0,Dy=0,Ds=0,Dxx=0,Dyy=0,Dss=0,Dxy=0,Dxs=0,Dys=0 ;
        float A [3*3], b [3] ;

        int dx = 0 ;
        int dy = 0 ;

        int iter, i, j ;

        for (iter = 0 ; iter < 5 ; ++iter) {

            x += dx ;
            y += dy ;

            pt = log
                + xo * x
                + yo * y
                + so * (s - s_min) ;


            /** @brief Index GSS @internal */
#define at(dx,dy,ds) (*( pt + (dx)*xo + (dy)*yo + (ds)*so))

            /** @brief Index matrix A @internal */
#define Aat(i,j)     (A[(i)+(j)*3])

            /* compute the gradient */
            Dx = 0.5 * (at(+1,0,0) - at(-1,0,0)) ;
            Dy = 0.5 * (at(0,+1,0) - at(0,-1,0));
            Ds = 0.5 * (at(0,0,+1) - at(0,0,-1)) ;

            /* compute the Hessian */
            Dxx = (at(+1,0,0) + at(-1,0,0) - 2.0 * at(0,0,0)) ;
            Dyy = (at(0,+1,0) + at(0,-1,0) - 2.0 * at(0,0,0)) ;
            Dss = (at(0,0,+1) + at(0,0,-1) - 2.0 * at(0,0,0)) ;

            Dxy = 0.25 * ( at(+1,+1,0) + at(-1,-1,0) - at(-1,+1,0) - at(+1,-1,0) ) ;
            Dxs = 0.25 * ( at(+1,0,+1) + at(-1,0,-1) - at(-1,0,+1) - at(+1,0,-1) ) ;
            Dys = 0.25 * ( at(0,+1,+1) + at(0,-1,-1) - at(0,-1,+1) - at(0,+1,-1) ) ;

            /* solve linear system ....................................... */
            Aat(0,0) = Dxx ;
            Aat(1,1) = Dyy ;
            Aat(2,2) = Dss ;
            Aat(0,1) = Aat(1,0) = Dxy ;
            Aat(0,2) = Aat(2,0) = Dxs ;
            Aat(1,2) = Aat(2,1) = Dys ;

            b[0] = - Dx ;
            b[1] = - Dy ;
            b[2] = - Ds ;

            /* Gauss elimination */
            for(j = 0 ; j < 3 ; ++j) {
                float maxa    = 0 ;
                float maxabsa = 0 ;
                int    maxi    = -1 ;
                float tmp ;

                /* look for the maximally stable pivot */
                for (i = j ; i < 3 ; ++i) {
                    float a    = Aat (i,j) ;
                    float absa = fabs (a) ;
                    if (absa > maxabsa) {
                        maxa    = a ;
                        maxabsa = absa ;
                        maxi    = i ;
                    }
                }

                /* if singular give up */
                if (maxabsa < 1e-10f) {
                    b[0] = 0 ;
                    b[1] = 0 ;
                    b[2] = 0 ;
                    break ;
                }

                i = maxi ;

                /* swap j-th row with i-th row and normalize j-th row */
                for(jj = j ; jj < 3 ; ++jj) {
                    tmp = Aat(i,jj) ; Aat(i,jj) = Aat(j,jj) ; Aat(j,jj) = tmp ;
                    Aat(j,jj) /= maxa ;
                }
                tmp = b[j] ; b[j] = b[i] ; b[i] = tmp ;
                b[j] /= maxa ;

                /* elimination */
                for (ii = j+1 ; ii < 3 ; ++ii) {
                    float x = Aat(ii,j) ;
                    for (jj = j ; jj < 3 ; ++jj) {
                        Aat(ii,jj) -= x * Aat(j,jj) ;
                    }
                    b[ii] -= x * b[j] ;
                }
            }

            /* backward substitution */
            for (i = 2 ; i > 0 ; --i) {
                float x = b[i] ;
                for (ii = i-1 ; ii >= 0 ; --ii) {
                    b[ii] -= x * Aat(ii,i) ;
                }
            }

            /* .......................................................... */
            /* If the translation of the keypoint is big, move the keypoint
            * and re-iterate the computation. Otherwise we are all set.
            */

            dx= ((b[0] >  0.6 && x < w - 2 - pad) ?  1 : 0)
                + ((b[0] < -0.6 && x > 1 + pad    ) ? -1 : 0) ;

            dy= ((b[1] >  0.6 && y < h - 2 - pad) ?  1 : 0)
                + ((b[1] < -0.6 && y > 1 + pad   ) ? -1 : 0) ;

            if (dx == 0 && dy == 0) break ;
        }

        /* check threshold and other conditions */
        {
            float val   = at(0,0,0)
                + 0.5 * (Dx * b[0] + Dy * b[1] + Ds * b[2]) ;
            float score = (Dxx+Dyy)*(Dxx+Dyy) / (Dxx*Dyy - Dxy*Dxy) ;
            float xn = x + b[0] ;
            float yn = y + b[1] ;
            float sn = s + b[2] ;

            int good =
                fabs (val)  > tp                  &&
                score           < (te+1)*(te+1)/te    &&
                score           >= 0                  &&
                fabs (b[0]) <  1.5                &&
                fabs (b[1]) <  1.5                &&
                fabs (b[2]) <  1.5                &&
                xn              >= 0                  &&
                xn              <= w - 1              &&
                yn              >= 0                  &&
                yn              <= h - 1              &&
                sn              >= s_min              &&
                sn              <= s_max ;


            if (good) {
                k-> o     = f->o_cur ;
                k-> ix    = x;
                k-> iy    = y;
                k-> is    = s;
                k-> s     = (float)sn ;
                k-> x     = (float)xn;
                k-> y     = (float)yn;
                k-> ix	  = (int)(k->x + 0.5);
                k-> iy	  = (int)(k->y + 0.5);
                k-> is    = (int)(k->s + 0.5); if(k->is < 0) k->is = 0; if(k->is > 2) k->is = 2;
                k-> sigma = (float)(f->sigma0 * pow (2.0f, sn/f->S));
                k-> peak = fabs(val) * FACTOR_LOG;
                k-> curvRatio = score;
#ifdef ADD_ATTRIBUTE
                k->log_dx = Dx;
                k->log_dy = Dy;
                k->log_ds = Ds;
                k->log_dxx = Dxx;
                k->log_dyy = Dyy;
                k->log_dss = Dss;
                k->log_dxy = Dxy;
                k->log_dxs = Dxs;
                k->log_dys = Dys;
                k->log_hessian_trace = Dxx+Dyy;
                k->log_hessian_det = Dxx*Dyy - Dxy*Dxy;
                k->log_hessian_eign1 = 0.5 * (k->log_hessian_trace + sqrt(k->log_hessian_trace * k->log_hessian_trace - 4 * k->log_hessian_det));
                k->log_hessian_eign2 = 0.5 * (k->log_hessian_trace - sqrt(k->log_hessian_trace * k->log_hessian_trace - 4 * k->log_hessian_det));
#endif
                ++ k ;
            }

        } /* done checking */
    } /* next keypoint to refine */

    /* update keypoint count */
    f-> nkeys = (int)(k - f->keys) ;
    return f->nkeys - st_key;
}
int keypoint_block_desc(SiftFilt *f, int blockid, int block_width, int block_height, FeatureList & featurelist, float xper)
{
    int cnt = 0;
    int i, iter, iori;
    SiftKeypoint * k;
    int block_row, block_col, cur_block_row, cur_block_col;
    int s_cur;
    gs_type * blockPtr = NULL;
    gs_type * p, *pt;
    float * src;
    int xi, yi, ys, xs, x1, x2, y1, y2, x0, y0;
    int left, right, up, bottom;
    float x, y, r2;
    int sig_x1, sig_y1, sig_x2, sig_y2, sig_x0, sig_y0;

    enum {nbins = 36} ;
    float hist [nbins], maxh ;
    float bin2Pi = (float)(nbins * 1.0 / CDVS_PI2);
    float const magnif = f-> magnif ;
    int nangles;
    float orientation[4];
    float rbuf [128] ;
    float hist_desc[NBP+2][NBP+2][NBO+2];

    float st0, ct0, SBP, SBP_1, sigma2;
    int W_Des, W_Ori;
    float wsigma = (float)(1.0 / (f->windowSize * f->windowSize * 2));

    block_row = blockid / f->block_row_n; //current block row
    block_col = blockid % f->block_row_n; //current block col
    //grad = (float *)vl_malloc(sizeof(float) * size * 2);

    Feature d;
    for(i = 0; i < f->nkeys; i++)
    {
        k = f->keys + i;
        if(k->mode == MODE_FINISH)
            continue;
        //the current block
        if(k->block_id == blockid){
            if((k->mode & MODE_RIGHT) || (k->mode & MODE_BOTTOM))
                continue;
        }
        //the left block
        else if(block_col > 0 && k->block_id + 1 == blockid ){
            if(k->mode & MODE_BOTTOM)
                continue;
        }
        //the upper block
        else if(block_row > 0 && k->block_id + f->block_row_n == blockid){
            if((k->mode & MODE_BOTTOM) && (k->mode & MODE_RIGHT))
                continue;
        }
        //the upper-left block
        else if(block_row > 0 && block_col > 0 && k->block_id + 1 + f->block_row_n == blockid){

        }
        else continue;

        cur_block_row = k->block_id / f->block_row_n; //current block row
        cur_block_col = k->block_id % f->block_row_n; //current block col
        s_cur = k->is;
        blockPtr = f->block_buffer + s_cur * f->block_buffer_len * BLOCK_SIZE;
        //////////////////////////////////////////////////////////////////////////
        SBP      = magnif * k->sigma + VL_EPSILON_D ;
        SBP_1	 = 1.0 / SBP;
        sigma2	 = 1.0 / (4.5 * k->sigma * k->sigma);
        W_Des    = (int)(floor(sqrt(2.0) * SBP * (NBP + 1) / 2.0 + 0.5));
        W_Ori	 = (int)(VL_MAX(floor(4.5 * k->sigma), 1));
        float W_Ori2 = W_Ori * W_Ori + 0.6f;

        k->mode = MODE_FINISH;
        memset(hist, 0, sizeof(hist));
        xi = (int) (k->x + 0.5);
        yi = (int) (k->y + 0.5);
        if(xi == BLOCK_WIDTH) xi = BLOCK_WIDTH - 1;	//risk
        if(yi == BLOCK_WIDTH) yi = BLOCK_WIDTH - 1;
        up = yi - W_Des;
        if(up <= 0 && cur_block_row == 0)
        {
            up = 1;
        }
        bottom = yi + W_Des;
        if(cur_block_row == f->block_col_n - 1 && bottom >= f->block_col_remain_height - 1)
        {
            bottom = f->block_col_remain_height - 2;
        }
        else if(cur_block_row == f->block_col_n - 2 && bottom >= f->block_col_remain_height - 1 + BLOCK_WIDTH)
        {
            bottom = f->block_col_remain_height - 2 + BLOCK_WIDTH;
        }
        left = xi - W_Des;
        if(left <= 0 && cur_block_col == 0)
        {
            left = 1;
        }
        right = xi + W_Des;
        if(cur_block_col == f->block_row_n - 1 && right >= f->block_row_remain_width - 1)
        {
            right = f->block_row_remain_width - 2;
        }
        else if(cur_block_col == f->block_row_n - 2 && right >= f->block_row_remain_width - 1 + BLOCK_WIDTH)
        {
            right = f->block_row_remain_width - 2 + BLOCK_WIDTH;
        }
        src = f->grad;
        for(ys = up; ys <= bottom; ys++)
        {
            y0 = ys;
            y1 = y0 - 1;
            y2 = y0 + 1;
            y = ys - k->y;
            y = y * y;

            if(y0 < 0){ sig_y0 = k->block_id-f->block_row_n; y0 += BLOCK_WIDTH;}
            else if(y0 < BLOCK_WIDTH){ sig_y0 = k->block_id;}
            else{ sig_y0 = k->block_id+f->block_row_n; y0 -= BLOCK_WIDTH;}

            if(y1 < 0){ sig_y1 = k->block_id-f->block_row_n; y1 += BLOCK_WIDTH;}
            else if(y1 < BLOCK_WIDTH){ sig_y1 = k->block_id;}
            else{ sig_y1 = k->block_id+f->block_row_n; y1 -= BLOCK_WIDTH;}

            if(y2 < 0){ sig_y2 = k->block_id-f->block_row_n; y2 += BLOCK_WIDTH;}
            else if(y2 < BLOCK_WIDTH){ sig_y2 = k->block_id;}
            else{ sig_y2 = k->block_id+f->block_row_n; y2 -= BLOCK_WIDTH;}

            y1 *= BLOCK_WIDTH;
            y2 *= BLOCK_WIDTH;
            y0 *= BLOCK_WIDTH;//*/

            for(xs = left; xs <= right; xs++)
            {
                float wgt, mod, ang, fbin, gx, gy, rbin;
                int bin, bin1, bin2;
                x0 = xs;
                x1 = x0 - 1;
                x2 = x0 + 1;
                x = xs - k->x;
                r2 = x * x + y;

                if(x0 < 0){ sig_x0 = -1; x0 += BLOCK_WIDTH;}
                else if(x0 < BLOCK_WIDTH){ sig_x0 = 0;}
                else{ sig_x0 = 1; x0 -= BLOCK_WIDTH;}

                if(x1 < 0){ sig_x1 = -1; x1 += BLOCK_WIDTH;}
                else if(x1 < BLOCK_WIDTH){ sig_x1 = 0;}
                else{ sig_x1 = 1; x1 -= BLOCK_WIDTH;}

                if(x2 < 0){ sig_x2 = -1; x2 += BLOCK_WIDTH;}
                else if(x2 < BLOCK_WIDTH){ sig_x2 = 0;}
                else{ sig_x2 = 1; x2 -= BLOCK_WIDTH;}//*/

                p = blockPtr + f->block_map[sig_y1 + sig_x0] + y1 + x0;
                pt = blockPtr + f->block_map[sig_y2 + sig_x0] + y2 + x0;
#ifdef ADD_ATTRIBUTE
                float x1y0, x2y0, x0y0, x0y1, x0y2, x1y1, x1y2, x2y1, x2y2;
                x0y1 = *p;
                x0y2 = *pt;
#endif
                gy = ((int)(*pt) - (int)(*p));

                p = blockPtr + f->block_map[sig_y0 + sig_x1] + y0 + x1;
                pt = blockPtr + f->block_map[sig_y0 + sig_x2] + y0 + x2;
                gx = ((int)(*pt) - (int)(*p));
#ifdef ADD_ATTRIBUTE
                if(xi == xs && yi == ys)
                {
                    x1y0 = *p;
                    x2y0 = *pt;
                    p = blockPtr + f->block_map[sig_y0 + sig_x0] + y0 + x0;
                    x0y0 = *p;
                    p = blockPtr + f->block_map[sig_y1 + sig_x1] + y1 + x1;
                    x1y1 = *p;
                    p = blockPtr + f->block_map[sig_y2 + sig_x1] + y2 + x1;
                    x1y2 = *p;
                    p = blockPtr + f->block_map[sig_y1 + sig_x2] + y1 + x2;
                    x2y1 = *p;
                    p = blockPtr + f->block_map[sig_y2 + sig_x2] + y2 + x2;
                    x2y2 = *p;
                    /* compute the gradient */
                    float Dx = 0.5 * x2y0 - x1y0;
                    float Dy = 0.5 * x0y2 - x0y1;
                    /* compute the Hessian */
                    float Dxx = x2y0 + x1y0 - 2.0 * x0y0 ;
                    float Dyy = x0y2 + x0y1 - 2.0 * x0y0 ;
                    float Dxy = 0.25 * ( x2y2 + x1y1 - x1y2 - x2y1 ) ;

                    k->g_dx = Dx;
                    k->g_dy = Dy;
                    k->g_dxx = Dxx;
                    k->g_dyy = Dyy;
                    k->g_dxy = Dxy;
                    k->g_hessian_trace = Dxx+Dyy;
                    k->g_hessian_det = Dxx*Dyy - Dxy*Dxy;
                    k->g_hessian_eign1 = 0.5 * (k->g_hessian_trace + sqrt(k->g_hessian_trace * k->g_hessian_trace - 4 * k->g_hessian_det));
                    k->g_hessian_eign2 = 0.5 * (k->g_hessian_trace - sqrt(k->g_hessian_trace * k->g_hessian_trace - 4 * k->g_hessian_det));
                }
#endif
                mod  = sqrt (gx*gx + gy*gy);cnt++;
                ang  = vl_fast_atan2_f (gy, gx);cnt++;
                *src++ = (float)mod;
                *src++ = (float)ang;

                if (r2 < W_Ori2)
                {
                    wgt  = fast_expn (r2 * sigma2) * mod;cnt++;
                    fbin = bin2Pi * ang - 0.5;
                    if(fbin < 0)
                    {
                        bin = -1;
                        bin1 = nbins - 1;
                        bin2 = 0;
                    }
                    else
                    {
                        bin = (int)fbin;
                        bin1 = bin;
                        bin2 = bin + 1;
                        if(bin2 >= nbins)
                        {
                            bin2 = 0;
                        }
                    }
                    rbin = (fbin - bin) * wgt ;
                    hist [bin1] += wgt - rbin ;
                    hist [bin2] += rbin ;
                }
            }//end x
        }//end y

        for (iter = 0; iter < 6; iter ++) {
            float prev  = hist [nbins - 1] ;
            float first = hist [0] ;
            int i ;
            for (i = 0; i < nbins - 1; i++) {
                float newh = (prev + hist[i] + hist[i+1]) / 3.0;
                prev = hist[i] ;
                hist[i] = newh ;
            }
            hist[i] = (prev + hist[i] + first) / 3.0 ;
        }

        maxh = 0 ;
        for (iori = 0 ; iori < nbins ; ++iori)
            maxh = VL_MAX (maxh, hist [iori]) ;

        /* find peaks within 80% from max */
        nangles = 0 ;
        for(iori = 0 ; iori < nbins ; ++iori) {
            float h0 = hist [iori] ;
            float hm = hist [(iori - 1 + nbins) % nbins] ;
            float hp = hist [(iori + 1 + nbins) % nbins] ;

            if (h0 > 0.8*maxh && h0 > hm && h0 > hp) {
                float di = - 0.5 * (hp - hm) / (hp + hm - 2 * h0) ;
                float th = CDVS_PI2 * (iori + di + 0.5) / nbins ;
                orientation [ nangles++ ] = th ;
                if( nangles == 4 )
                    break;
            }
        }//*/
        //////////////////////////////////////////////////////////////////////////
        //output the detector's character
        d.x = (k->x + BLOCK_WIDTH * cur_block_col) * xper;
        d.y = (k->y + BLOCK_WIDTH * cur_block_row) * xper;
        d.scale = k->sigma * xper;
        d.peak = k->peak;
        d.iscale = k->is;
        d.octave = k->o;
        d.curvRatio = k->curvRatio;
#ifdef ADD_ATTRIBUTE
        d.g_dx = k->g_dx;
        d.g_dy = k->g_dy;
        d.g_dxx = k->g_dxx;
        d.g_dyy = k->g_dyy;
        d.g_dxy = k->g_dxy;
        d.g_hessian_det = k->g_hessian_det;
        d.g_hessian_eign1 = k->g_hessian_eign1;
        d.g_hessian_eign2 = k->g_hessian_eign2;
        d.g_hessian_trace = k->g_hessian_trace;
        d.log_dx = k->log_dx;
        d.log_dy = k->log_dy;
        d.log_ds = k->log_ds;
        d.log_dxx = k->log_dxx;
        d.log_dyy = k->log_dyy;
        d.log_dss = k->log_dss;
        d.log_dxy = k->log_dxy;
        d.log_dxs = k->log_dxs;
        d.log_dys = k->log_dys;
        d.log_hessian_det = k->log_hessian_det;
        d.log_hessian_eign1 = k->log_hessian_eign1;
        d.log_hessian_eign2 = k->log_hessian_eign2;
        d.log_hessian_trace = k->log_hessian_trace;
#endif
        //////////////////////////////////////////////////////////////////////////
        int binx, biny, bint, binx1, biny1, bint1;
        for(iori = 0; iori < nangles; iori++)
        {
            int bin;
            float *dpt ;
            float norm;
            float angle0 = orientation[iori];
            st0 = sin (angle0) * SBP_1;cnt++;
            ct0 = cos (angle0) * SBP_1;cnt++;
            d.orientation = angle0;
            memset(hist_desc, 0, sizeof(hist_desc));
            src = f->grad;
            float dx, dy, nx, ny, nx1, ny1, nx2, ny2;
            float rbinx, rbiny, rbint;
            float v_r1, v_r0, v_rc11, v_rc10, v_rc01, v_rc00;
            float v_rco111, v_rco110, v_rco101, v_rco100;
            float v_rco011, v_rco010, v_rco001, v_rco000;
            dy = up - k->y - 1;
            dx = left - k->x - 1;
            nx2 = ct0 * dx + st0 * dy;
            ny2 = -st0 * dx + ct0 * dy;
            for(y0 = up; y0 <= bottom; y0++)
            {
                //float dy = y0 - k->y;
                nx1 = (nx2 += st0);
                ny1 = (ny2 += ct0);
                for(x0 = left; x0 <= right; x0++)
                {
                    nx1 += ct0;
                    ny1 -= st0;
                    if(nx1 >= 2.5 || nx1 < -2.5 || ny1 >= 2.5 || ny1 < -2.5) { src += 2; continue;}
                    //get the gradient value
                    float mod = *src++;
                    float angle = *src++;

                    float theta = angle - angle0;
                    if(theta >= CDVS_PI2) theta -= CDVS_PI2;
                    if(theta < 0) theta += CDVS_PI2;

                    float nt = BPO_P_PI2 * theta;

                    //get the weight
                    //float win = fast_expn((nx*nx + ny*ny) * wsigma) * mod;
                    float win = fast_expn((nx1*nx1 + ny1*ny1) * wsigma) * mod;cnt++;

                    //apply fast tri-linear interpolation
                    //nx += 2.5;
                    //ny += 2.5;
                    nx = nx1 + 2.5;
                    ny = ny1 + 2.5;
                    //binx = vl_floor_f (nx) ;
                    //biny = vl_floor_f (ny) ;
                    //bint = vl_floor_f (nt) ;
                    binx = (int)nx; biny = (int)ny; bint = (int)nt;
                    binx1 = binx + 1;
                    biny1 = biny + 1;
                    bint1 = bint + 1;
                    rbinx = nx - binx ;
                    rbiny = ny - biny ;
                    rbint = nt - bint ;
                    /*int         dbinx ;
                    int         dbiny ;
                    int         dbint ;//*/

                    v_r1 = win * rbinx, v_r0 = win - v_r1;
                    v_rc11 = v_r1 * rbiny, v_rc10 = v_r1 - v_rc11;
                    v_rc01 = v_r0 * rbiny, v_rc00 = v_r0 - v_rc01;
                    v_rco111 = v_rc11 * rbint, v_rco110 = v_rc11 - v_rco111;
                    v_rco101 = v_rc10 * rbint, v_rco100 = v_rc10 - v_rco101;
                    v_rco011 = v_rc01 * rbint, v_rco010 = v_rc01 - v_rco011;
                    v_rco001 = v_rc00 * rbint, v_rco000 = v_rc00 - v_rco001;

                    hist_desc[binx][biny][bint] += v_rco000;
                    hist_desc[binx][biny][bint1] += v_rco001;
                    hist_desc[binx][biny1][bint] += v_rco010;
                    hist_desc[binx][biny1][bint1] += v_rco011;
                    hist_desc[binx1][biny][bint] += v_rco100;
                    hist_desc[binx1][biny][bint1] += v_rco101;
                    hist_desc[binx1][biny1][bint] += v_rco110;
                    hist_desc[binx1][biny1][bint1] += v_rco111;  //*/

                }//end x
            }//end y

            //output the descriptor
            dpt = rbuf;
            for(biny = 1; biny <= NBP; biny++)
            {
                for(binx = 1; binx <= NBP; binx++)
                {
                    hist_desc[binx][biny][0] += hist_desc[binx][biny][NBO];
                    hist_desc[binx][biny][1] += hist_desc[binx][biny][NBO+1];
                    for(bint = 0; bint < NBO; bint++)
                    {
                        *dpt++ = hist_desc[binx][biny][bint];
                    }
                }
            }//*/
            /* Normalize the histogram to L2 unit length. */
            norm = normalize_histogram (rbuf) ;

            /* Set the descriptor to zero if it is lower than our norm_threshold */
            if(f-> norm_thresh && norm < f-> norm_thresh) {
                for(bin = 0; bin < 128; ++ bin)
                    rbuf[bin] = 0;
            }
            else {

                /* Truncate at 0.2. */
                for(bin = 0; bin < 128; ++ bin) {
                    if (rbuf[bin] > 0.2) rbuf[bin] = 0.2;
                }

                /* Normalize again. */
                normalize_histogram (rbuf) ;
            }
            for (int j = 0 ; j < 128 ; ++j) {
                float x = 512.0F * rbuf [j] ;
                x = (x < 255.0F) ? x : 255.0F ;
                d.descr[j] = (unsigned char) x ;		// return uint8
            }
            featurelist.append(d);
        }//end orient
    }//end key
    return cnt;
}
