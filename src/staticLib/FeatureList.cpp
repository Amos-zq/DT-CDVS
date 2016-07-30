#include "defines.h"
#include "FeatureList.h"
#include "CsscCoordinateCoding.h"
#include "LookupTable.h"
extern LookupTable lookupTable;

//for create featurelist from image 
FeatureList::FeatureList(int targetLength):features(),numExtracted(0),originalMaxResolution(0),imageHeight(0),imageWidth(0) {
    if(targetLength >= 16384)       qdescr_size = 128;
    else if(targetLength >= 8192)   qdescr_size = 80;
    else if(targetLength >= 4096)   qdescr_size = 64;
    else if(targetLength >= 2048)   qdescr_size = 40;
    else                            qdescr_size = 20;
}
void FeatureList::copyFeatureList(FeatureList& featureList) {
    imageWidth = featureList.imageWidth;
    imageHeight = featureList.imageHeight;
    originalMaxResolution = featureList.originalMaxResolution;
    qdescr_size = featureList.qdescr_size;
    numExtracted = featureList.numExtracted;
    for(int i = 0; i < numExtracted; ++i) {
        features.push_back(featureList.features[i]);
    }
}
FeatureList::~FeatureList() {
    std::vector<Feature> featureNull = std::vector<Feature> ();
    features.swap(featureNull);
}

void FeatureList::setDescLength(size_t descLength) {
    if(descLength >= 16384)       qdescr_size = 128;
    else if(descLength >= 8192)   qdescr_size = 80;
    else if(descLength >= 4096)   qdescr_size = 64;
    else if(descLength >= 2048)   qdescr_size = 40;
    else                          qdescr_size = 20;
}
void FeatureList::setResolution(int imgWidth, int imgHeight, int originalMaxRes) {
    imageHeight = imgHeight;
    imageWidth = imgWidth;
    originalMaxResolution = originalMaxRes;
}
void FeatureList::setRelevantPoints(int num) {
    if (num == 0)   return;
    int numNotRelevantPoints = CDVS_MAX(int(numExtracted) - num, 0);

    for(std::vector<Feature>::iterator p=features.begin(); p<features.end()-numNotRelevantPoints; ++p) {
        p->relevance = 1;
    }
}

bool FeatureList::sortSpatialIndex(const Feature &f1, const Feature &f2) {
    return f1.spatialIndex < f2.spatialIndex;
}
bool FeatureList::sortRelevance(const Feature &f1, const Feature &f2) {
    return f1.relevance > f2.relevance;
}
void FeatureList::sortFeatures(bool (*compare)(const Feature&, const Feature&)) {
    stable_sort(features.begin(), features.end(), compare);
}

size_t FeatureList::NumExtracted() const {
    return numExtracted;
}
size_t FeatureList::size() {
    size_t ret = 0;
    ret += sizeof(imageWidth);
    ret += sizeof(imageHeight);
    ret += sizeof(originalMaxResolution);
    ret += sizeof(unsigned int);
    ret += (sizeof(float) * 136 + sizeof(int) + sizeof(unsigned short)) * numExtracted;
    return ret;
}

void FeatureList::append(const Feature& feature) {
    features.push_back(feature);
    ++numExtracted;
}
void FeatureList::selectFirst(int n) {
    if (n < numExtracted) {
        features.erase(features.begin() + n, features.end());
        numExtracted = n;
    }
}

void FeatureList::readSiftStream(BitInputStream &in) {
    in.read((unsigned char*)&imageWidth, sizeof(imageWidth) * 8);
    in.read((unsigned char*)&imageHeight, sizeof(imageHeight) * 8);
    in.read((unsigned char*)&originalMaxResolution, sizeof(originalMaxResolution) * 8);
    in.read((unsigned char*)&numExtracted, sizeof(unsigned int) * 8);
    features.resize(numExtracted, Feature());
    for(unsigned int i=0; i<numExtracted; i++) {
        features[i].fromStream(in);
    }
}
void FeatureList::writeSiftStream(BitOutputStream &out) const{
    out.write((unsigned char*)&imageWidth, sizeof(imageWidth) * 8);
    out.write((unsigned char*)&imageHeight, sizeof(imageHeight) * 8);
    out.write((unsigned char*)&originalMaxResolution, sizeof(originalMaxResolution) * 8);
    out.write((unsigned char*)&numExtracted, sizeof(unsigned int) * 8);
    for(unsigned int i=0; i<numExtracted; i++) {
        features[i].toStream(out);
    }
}
void FeatureList::writeSift(const char *filename) const {
    FILE * file = fopen(filename, "wb");
    fwrite(&imageWidth, sizeof(imageWidth), 1, file);
    fwrite(&imageHeight, sizeof(imageHeight), 1, file);
    fwrite(&originalMaxResolution, sizeof(originalMaxResolution), 1, file);
    fwrite(&numExtracted, sizeof(unsigned int), 1, file);
    for(unsigned int i=0; i<numExtracted; i++) {
        features[i].toFile(file);
    }
    fclose(file);
}
void FeatureList::readSift(const char *filename) {
    FILE * file = fopen(filename, "rb");
    fread(&imageWidth, sizeof(imageWidth), 1, file);
    fread(&imageHeight, sizeof(imageHeight), 1, file);
    fread(&originalMaxResolution, sizeof(originalMaxResolution), 1, file);
    fread(&numExtracted, sizeof(unsigned int), 1, file);
    features.resize(numExtracted, Feature());
    for(unsigned int i=0; i<numExtracted; i++) {
        features[i].fromFile(file);
    }
    fclose(file);
}

void FeatureList::compressHM(int targetLength) {
    int N = int(qdescr_size >> 2);
    for(std::vector<Feature>::iterator f=features.begin(); f<features.end(); ++f) {
        if (f->descr[0] == -1.0f)   throw CdvsException("FeatureList::compressHM error: data not valid");
        int i = 0,  j = 0;
        for(i = 0; i < N; i++ ) {
            int group_id    = LookupTable::priority_list[i][0];
            int element_id  = LookupTable::priority_list[i][1];
            int v0 = 0, v1 = 0, v2 = 0, v3 = 0;

            int hist_pos0 = LookupTable::histogram_groups[group_id][0] << 3;
            int hist_pos1 = LookupTable::histogram_groups[group_id][1] << 3;
            int hist_pos2 = LookupTable::histogram_groups[group_id][2] << 3;
            int hist_pos3 = LookupTable::histogram_groups[group_id][3] << 3;

            float *d0 = f->descr + hist_pos0;
            float *d1 = f->descr + hist_pos1;
            float *d2 = f->descr + hist_pos2;
            float *d3 = f->descr + hist_pos3;

            switch(element_id) {
            case 0:
                // transform A
                v0 = (int)(d0[2] - d0[6])/2;    v1 = (int)(d1[2] - d1[6])/2;
                // transform B
                v2 = (int)(d2[0] - d2[4])/2;    v3 = (int)(d3[0] - d3[4])/2;
                break;
            case 1:
                // transform A
                v0 = (int)(d0[3] - d0[7])/2;    v1 = (int)(d1[3] - d1[7])/2;
                // transform B
                v2 = (int)(d2[1] - d2[5])/2;    v3 = (int)(d3[1] - d3[5])/2;
                break;
            case 2:
                // transform A
                v0 = (int)(d0[0] - d0[1])/2;    v1 = (int)(d1[0] - d1[1])/2;
                // transform B
                v2 = (int)(d2[7] - d2[0])/2;    v3 = (int)(d3[7] - d3[0])/2;
                break;
            case 3:
                // transform A
                v0 = (int)(d0[2] - d0[3])/2;    v1 = (int)(d1[2] - d1[3])/2;
                // transform B
                v2 = (int)(d2[1] - d2[2])/2;    v3 = (int)(d3[1] - d3[2])/2;
                break;
            case 4:
                // transform A
                v0 = (int)(d0[4] - d0[5])/2;    v1 = (int)(d1[4] - d1[5])/2;
                // transform B
                v2 = (int)(d2[3] - d2[4])/2;    v3 = (int)(d3[3] - d3[4])/2;
                break;
            case 5:
                // transform A
                v0 = (int)(d0[6] - d0[7])/2;    v1 = (int)(d1[6] - d1[7])/2;
                // transform B
                v2 = (int)(d2[5] - d2[6])/2;    v3 = (int)(d3[5] - d3[6])/2;
                break;
            case 6:
                // transform A
                v0 = (int)(d0[0] + d0[4] - d0[2] - d0[6])/4;    v1 = (int)(d1[0] + d1[4] - d1[2] - d1[6])/4;
                // transform B
                v2 = (int)(d2[1] + d2[5] - d2[3] - d2[7])/4;    v3 = (int)(d3[1] + d3[5] - d3[3] - d3[7])/4;
                break;
            case 7:
                // transform A
                v0 = (int)(d0[0] + d0[2]  + d0[4] + d0[6] - d0[1] - d0[3] - d0[5] - d0[7])/8;
                v1 = (int)(d1[0] + d1[2]  + d1[4] + d1[6] - d1[1] - d1[3] - d1[5] - d1[7])/8;
                // transform B
                v2 = (int)(d2[0] + d2[1]  + d2[2] + d2[3] - d2[4] - d2[5] - d2[6] - d2[7])/8;
                v3 = (int)(d3[0] + d3[1]  + d3[2] + d3[3] - d3[4] - d3[5] - d3[6] - d3[7])/8;
                break;
            }

            if(v0 > LookupTable::hm_ter_quant_table[hist_pos0 + element_id][1])         f->qdescr[j++] = 2;
            else if(v0 > LookupTable::hm_ter_quant_table[hist_pos0 + element_id][0])    f->qdescr[j++] = 1;
            else                                                                        f->qdescr[j++] = 0;

            if(v1 > LookupTable::hm_ter_quant_table[hist_pos1 + element_id][1])         f->qdescr[j++] = 2;
            else if(v1 > LookupTable::hm_ter_quant_table[hist_pos1 + element_id][0])    f->qdescr[j++] = 1;
            else                                                                        f->qdescr[j++] = 0;

            if(v2 > LookupTable::hm_ter_quant_table[hist_pos2 + element_id][1])         f->qdescr[j++] = 2;
            else if(v2 > LookupTable::hm_ter_quant_table[hist_pos2 + element_id][0])    f->qdescr[j++] = 1;
            else                                                                        f->qdescr[j++] = 0;

            if(v3 > LookupTable::hm_ter_quant_table[hist_pos3 + element_id][1])         f->qdescr[j++] = 2;
            else if(v3 > LookupTable::hm_ter_quant_table[hist_pos3 + element_id][0])    f->qdescr[j++] = 1;
            else                                                                        f->qdescr[j++] = 0;
        }
    }
}
void FeatureList::writeStream(BitOutputStream &writer, bool writeRelevance, int nFeatures) {
    if(numExtracted <= 0)     return;
    if(nFeatures > numExtracted)    nFeatures = int(numExtracted);
    writer.write(qdescr_size >> 2, 6);
#ifdef USE_DESCR_ARITHMETIC_CODER
    AC_encoder ace;
    AC_model acm;
    ace.init(writer);
    acm.init(AC_MAX_VALUE4, 0, 1);
    for(std::vector<Feature>::iterator f=features.begin(); f < features.begin() + nFeatures; ++f) {
        int vals[4];
        int val_ind = 0;
        for(int i=0; i<qdescr_size; i++) {
            vals[val_ind] = f->qdescr[i];
            val_ind++;
            if(val_ind == 4) {
                unsigned int ac_val = vals[0] + 3*vals[1] + 9*vals[2] + 27*vals[3];
                ace.encode_symbol(acm, ac_val);
                vals[0] = vals[1] = vals[2] = vals[3] = 0;
                val_ind = 0;
            }
        }
    }
    ace.done();
    acm.done();
#else
    for(std::vector<Feature>::iterator f=features.begin(); f < features.begin() + nFeatures; ++f) {
        for(unsigned int i=0; i < qdescr_size; i++) {
            if(f->qdescr[i] == 1) {
                writer.write(0);// 1 => '0'
            } else {
                writer.write(1);
                if(f->qdescr[i] == 0) {
                    writer.write(0);// 0 => '10'
                } else {
                    writer.write(1);// 2 => '11'
                }
            }
        }
    }
#endif // USE_DESCR_ARITHMETIC_CODER

    if(writeRelevance)  {
        for(std::vector<Feature>::iterator f=features.begin(); f<features.begin() + nFeatures; ++f)
            writer.write(f->relevance>0);
    }
}
void FeatureList::readStream(BitInputStream &reader, bool readRelevance) {
    if(numExtracted == 0)   return;
    qdescr_size = reader.read(6) * 4;
    if(qdescr_size != 20 && qdescr_size != 40 && qdescr_size != 64 && qdescr_size != 80 && qdescr_size != 128) {
        qdescr_size = 0;
        return;
    }
#ifdef USE_DESCR_ARITHMETIC_CODER
    AC_decoder acd;
    AC_model acm;
    acd.init(reader);
    acm.init(AC_MAX_VALUE4, 0, 1);

    for(std::vector<Feature>::iterator f=features.begin(); f<features.end(); ++f) {
        unsigned char uc;
        int read_elems = qdescr_size;
        int byte_offset = 0, descr_bytes = 0;
        int vals[4];

        int i = 0;
        while(read_elems) {
            unsigned int ac_val;
            ac_val = acd.decode_symbol(acm);

            vals[3] = ac_val/27;
            ac_val -= vals[3]*27;
            vals[2] = ac_val/9;
            ac_val -= vals[2]*9;
            vals[1] = ac_val/3;
            ac_val -= vals[1]*3;
            vals[0] = ac_val;
            f->qdescr[i++] = vals[0];
            f->qdescr[i++] = vals[1];
            f->qdescr[i++] = vals[2];
            f->qdescr[i++] = vals[3];
            read_elems-=4;
            if(read_elems == 0) break;
        }
        for(i = 0; i < 128; i++) {
            f->descr[i] = 0;
        }

        // Pack descriptor elements into bytes: 4 elements per byte.
        for(i = 0; i < qdescr_size; i++) {
            int val;
            if(byte_offset == 0) {
                f->ucdescr[descr_bytes] = 0;
                uc = 0;
            }
            val = f->qdescr[i];
            if(val == 2) val = 3;
            uc |= (val & 3) << byte_offset;
            byte_offset+=2;
            if(byte_offset == 8) {
                f->ucdescr[descr_bytes] = uc;
                byte_offset = 0;
                // go to next byte;
                descr_bytes++;
            }
        }
    }
    acd.done();
    acm.done();
#else
    for(std::vector<Feature>::iterator f=features.begin(); f<features.end(); ++f) {
        unsigned char uc = 0;
        int read_elems = qdescr_size;
        unsigned int i = 0;
        while(read_elems) {
            int v = reader.read();
            if(v == 0) {
                f->qdescr[i] = 1;
            } else {
                v = reader.read();
                if(v == 0) {
                    f->qdescr[i] = 0;
                } else {
                    f->qdescr[i] = 2;
                }
            }
            i++;
            read_elems--;
        }
        for(i = 0; i < 128; i++) f->descr[i] = 0;

        // Pack descriptor elements into bytes: 4 elements per byte.
        int byte_offset = 0, descr_bytes = 0;
        for(i = 0; i < qdescr_size; i++) {
            if(byte_offset == 0) {
                f->ucdescr[descr_bytes] = 0;
                uc = 0;
            }
            int val = f->qdescr[i];
            if(val == 2) val = 3;
            uc |= (val & 3) << byte_offset;
            byte_offset+=2;
            if(byte_offset == 8) {
                f->ucdescr[descr_bytes] = uc;
                byte_offset = 0;
                descr_bytes++;
            }
        }
    }
#endif // USE_DESCR_ARITHMETIC_CODER

    if(readRelevance) {
        for(std::vector<Feature>::iterator f=features.begin(); f<features.end(); ++f) {
            f->relevance = reader.read(1);
        }
    }
}

int FeatureList::matchSift(PointPairs &pairs, FeatureList &otherFeatureList, float ratioThreshold) {
    pairs.nMatched = 0;
    if(features.size()>1 && otherFeatureList.features.size()>1) {
        ratioThreshold = ratioThreshold*ratioThreshold;
        float distance, minDistance, secondMinDistance;
        int minDistanceInd, featureInd = 0, otherFeatureInd;
        Match match;
        std::vector<Match> matches;

        for(vector<Feature>::const_iterator f = features.begin(); f < features.end(); ++f) {
            // Find the two nearest descriptors contained in otherFeatureList and the relative distances
            minDistance = 1e10;
            secondMinDistance = 1e10;
            otherFeatureInd = 0;
            minDistanceInd = 0;

            for(vector<Feature>::const_iterator o = otherFeatureList.features.begin(); o < otherFeatureList.features.end(); ++o) {
                distance = 0;
                for(int i=0;i<128;++i) {
                    distance += (f->descr[i] - o->descr[i]) * (f->descr[i] - o->descr[i]);
                    if(distance >= secondMinDistance)   break;
                }

                if(distance < minDistance) {
                    secondMinDistance = minDistance;
                    minDistance = distance;
                    minDistanceInd = otherFeatureInd;
                } else {
                    if(distance<secondMinDistance) {
                        secondMinDistance = distance;
                    }
                }
                ++otherFeatureInd;
            }

            // If the ratio test is passed the indices of the features are saved
            if(secondMinDistance > 0 && minDistance <= ratioThreshold * secondMinDistance) {
                match.featureInd = featureInd;
                match.otherFeatureInd = minDistanceInd;
                match.weight = cos(HALF_PI * sqrt((double)minDistance/(double)secondMinDistance));
                matches.push_back(match);
            }
            ++featureInd;
        }

        //bi-unique matchings deleting possible repetitions in otherFeatureInd
        if(matches.size()>0) {
            sort(matches.begin(), matches.end(), Match::sortMatchByWeight);
            int lastFeatureIndex = -1;

            for(vector<Match>::const_iterator m = matches.begin(); m < matches.end(); ++m) {
                if(lastFeatureIndex != m->otherFeatureInd) {
                    lastFeatureIndex = m->otherFeatureInd;
                    pairs.addPair(features[m->featureInd], otherFeatureList.features[m->otherFeatureInd], m->featureInd, m->otherFeatureInd, m->weight);
                }
            }
        }
    }
    return pairs.nMatched;
}
int FeatureList::matchOneWay(PointPairs &pairs, FeatureList &otherFeatureList, float ratioThreshold) {
    pairs.nMatched = 0;
    if(features.size() > 1 && otherFeatureList.features.size() > 1) {
        ratioThreshold = ratioThreshold*ratioThreshold;
        int distance, minDistance, secondMinDistance;
        int minDistanceInd, featureInd = 0, otherFeatureInd;
        int matched_bytes = qdescr_size >> 2;
        if((otherFeatureList.qdescr_size >> 2) < matched_bytes) {
            matched_bytes = otherFeatureList.qdescr_size >> 2;
        }

        Match match;
        std::vector<Match> matches;

        for(vector<Feature>::const_iterator f = features.begin(); f < features.end(); ++f) {
            // Find the two neares descriptors contained in otherFeatureList and the relative distances
            minDistance = 65536;
            secondMinDistance = 65536;
            otherFeatureInd = 0;
            minDistanceInd = 0;
            unsigned char *ucdescr = (unsigned char *)f->ucdescr;

            for(vector<Feature>::const_iterator o = otherFeatureList.features.begin(); o < otherFeatureList.features.end(); ++o) {
                unsigned char *otherUcdescr = (unsigned char *)o->ucdescr;
                distance = 0;
                for(unsigned int i = 0;i < matched_bytes;++i) {
                    distance += lookupTable.fHanming[ucdescr[i] ^ otherUcdescr[i]];
                }

                if(distance<minDistance) {
                    secondMinDistance = minDistance;
                    minDistance = distance;
                    minDistanceInd = otherFeatureInd;
                } else {
                    if(distance<secondMinDistance) {
                        secondMinDistance = distance;
                    }
                }
                ++otherFeatureInd;
            }

            // If the ratio test is passed the indices of the features are saved
            if((float)minDistance <= ratioThreshold*(float)secondMinDistance && secondMinDistance>0) {
                match.featureInd = featureInd;
                match.otherFeatureInd = minDistanceInd;
                match.weight = cos(HALF_PI * sqrt((double)minDistance/(double)secondMinDistance));
                matches.push_back(match);
            }
            ++featureInd;
        }

        // bi-unique matchings deleting possible repetitions in otherFeatureInd
        if(matches.size()>0) {
            std::sort(matches.begin(), matches.end(), Match::sortMatchByWeight);
            int lastFeatureIndex = -1;
            for(std::vector<Match>::const_iterator m=matches.begin(); m<matches.end(); ++m) {
                if(lastFeatureIndex != m->otherFeatureInd) {
                    lastFeatureIndex = m->otherFeatureInd;
                    pairs.addPair(features[m->featureInd], otherFeatureList.features[m->otherFeatureInd], m->featureInd, m->otherFeatureInd, m->weight);
                }
            }
        }
    }
    return pairs.nMatched;
}
int FeatureList::matchTwoWay(PointPairs &pairs, FeatureList &otherFeatureList, float ratioThreshold) {
    pairs.nMatched = 0;
    if(features.size()>1 && otherFeatureList.features.size()>1)  {
        Match match;
        std::vector<Match> matches1;
        std::vector<Match> matches1clean;
        std::vector<Match> matches2;
        std::vector<Match> matches2clean;
        
        ratioThreshold *= ratioThreshold;
        int distance, minDistance, secondMinDistance;
        int minDistanceInd, featureInd = 0, otherFeatureInd;

        int** distanceMatrix = new int*[features.size()];
        for (unsigned int i=0; i < features.size(); i++)
            distanceMatrix[i] = new int[otherFeatureList.features.size()];


        int matched_bytes = qdescr_size >> 2;
        if((otherFeatureList.qdescr_size >> 2) < matched_bytes) {
            matched_bytes = otherFeatureList.qdescr_size >> 2;
        }

        for(vector<Feature>::const_iterator f = features.begin(); f < features.end(); ++f) {
            // Find the two neares descriptors contained in otherFeatureList and the relative distances
            minDistance = 65536;
            secondMinDistance = 65536;
            otherFeatureInd = 0;
            minDistanceInd = 0;
            unsigned char* ucdescr = (unsigned char *)f->ucdescr;

            for(vector<Feature>::const_iterator o = otherFeatureList.features.begin(); o < otherFeatureList.features.end(); ++o) {
                unsigned char *otherUcdescr = (unsigned char *)o->ucdescr;
                distance = 0;
                for(unsigned int i=0;i<matched_bytes;++i) {
                    distance += lookupTable.fHanming[ucdescr[i] ^ otherUcdescr[i]];
                }

                distanceMatrix[featureInd][otherFeatureInd] = distance;

                if(distance<minDistance) {
                    secondMinDistance = minDistance;
                    minDistance = distance;
                    minDistanceInd = otherFeatureInd;
                } else {
                    if(distance<secondMinDistance) {
                        secondMinDistance = distance;
                    }
                }
                ++otherFeatureInd;
            }

            // If the ratio test is passed the indices of the features are saved
            if((float)minDistance <= ratioThreshold * (float)secondMinDistance && secondMinDistance > 0) {
                match.featureInd = featureInd;
                match.otherFeatureInd = minDistanceInd;
                match.weight = (float)cos(HALF_PI * sqrt((float)minDistance / (float)secondMinDistance));
                matches1.push_back(match);
            }
            ++featureInd;
        }

        //bi-unique matchings deleting possible repetitions in otherFeatureInd
        if(matches1.size()>0) {
            std::sort(matches1.begin(), matches1.end(), Match::sortMatchByWeight);
            int lastFeatureIndex = -1;

            for(std::vector<Match>::const_iterator m=matches1.begin(); m<matches1.end(); ++m) {
                if(lastFeatureIndex != m->otherFeatureInd) {
                    lastFeatureIndex = m->otherFeatureInd;
                    Match tempMatch;
                    tempMatch.featureInd = m->featureInd;
                    tempMatch.otherFeatureInd = m->otherFeatureInd;
                    tempMatch.weight = m->weight;
                    matches1clean.push_back(tempMatch);
                }
            }
        }
        matches1.clear();

        // Select the two nearest descriptors
        otherFeatureInd = 0;
        for(vector<Feature>::const_iterator o = otherFeatureList.features.begin(); o < otherFeatureList.features.end(); ++o) {
            // Find the two nearest descriptors contained in otherFeatureList and the relative distances between f
            minDistance = 1e10;
            secondMinDistance = 1e10;
            featureInd = 0;
            minDistanceInd = 0;
            const unsigned char *ucdescr = o->ucdescr;
            for(vector<Feature>::const_iterator f = features.begin(); f < features.end(); ++f) {
                distance = distanceMatrix[featureInd][otherFeatureInd];
                if(distance<minDistance) {
                    secondMinDistance = minDistance;
                    minDistance = distance;
                    minDistanceInd = featureInd;
                } else {
                    if(distance<secondMinDistance) {
                        secondMinDistance = distance;
                    }
                }
                ++featureInd;
            }

            // If the ratio test is passed the indices of the features are saved
            if ((float)minDistance < ratioThreshold * (float)secondMinDistance && secondMinDistance>0) {
                match.featureInd = minDistanceInd;
                match.otherFeatureInd = otherFeatureInd;
                match.weight = (float)cos(HALF_PI * sqrt((float)minDistance / (float)secondMinDistance));
                matches2.push_back(match);
            }
            ++otherFeatureInd;
        }

        //bi-unique matchings deleting possible repetitions in otherFeatureInd
        if(matches2.size()>0) {
            std::sort(matches2.begin(), matches2.end(), Match::sortMatchByWeight);
            int lastFeatureIndex = -1;

            for(std::vector<Match>::const_iterator m=matches2.begin(); m<matches2.end(); ++m) {
                if(lastFeatureIndex != m->featureInd) {
                    lastFeatureIndex = m->featureInd;
                    Match tempMatch;
                    tempMatch.featureInd = m->featureInd;
                    tempMatch.otherFeatureInd = m->otherFeatureInd;
                    tempMatch.weight = m->weight;
                    matches2clean.push_back(tempMatch);
                }
            }
        }
        matches2.clear();

        for(std::vector<Match>::iterator m1 = matches1clean.begin(); m1 < matches1clean.end(); ++m1) {
            int flag=0;
            for(std::vector<Match>::iterator m2 = matches2clean.begin(); m2 < matches2clean.end(); ++m2) {
                if(m1->featureInd == m2->featureInd && m1->otherFeatureInd == m2->otherFeatureInd) {
                    pairs.addPair(features[m1->featureInd], otherFeatureList.features[m2->otherFeatureInd], m1->featureInd, m2->otherFeatureInd, m1->weight);
                    flag=1; break;
                }
            }
            if(flag==0) {
                pairs.addPair(features[m1->featureInd], otherFeatureList.features[m1->otherFeatureInd], m1->featureInd, m1->otherFeatureInd, m1->weight*0.5);
            }
        }

        for(std::vector<Match>::iterator m2 = matches2clean.begin(); m2 < matches2clean.end(); ++m2) {
            int flag=0;
            for(std::vector<Match>::iterator m1 = matches1clean.begin(); m1 < matches1clean.end(); ++m1) {
                if(m1->featureInd == m2->featureInd && m1->otherFeatureInd == m2->otherFeatureInd){
                    flag=1; break;
                }
            }
            if(flag==0) {
                pairs.addPair(features[m2->featureInd], otherFeatureList.features[m2->otherFeatureInd], m2->featureInd, m2->otherFeatureInd, m2->weight*0.5);
            }
        }

        for (unsigned int i=0; i<features.size(); i++)
            if (distanceMatrix[i] != NULL) delete [] distanceMatrix[i];
        delete [] distanceMatrix;
    }
    return pairs.nMatched;

}

int FeatureList::computeMaxPoints(const Parameters & params, int targetBits) {
    static const double bitsPerKeypoint[] = {7.5, 7.0, 5.0, 4.8, 4.7, 4.6, 4.5};
    // compute a rough estimate keypoint number (num1)
    double availableBits = targetBits;				// compute the available bits
    double estimate_bits = qdescr_size * 1.66666 + bitsPerKeypoint[params.modeId - 7];	// bits per key point to encode location bits
    int num1 = (int)(availableBits / estimate_bits + 0.5);	// first approximation of the number of keypoints
    if (num1 > numExtracted)
        num1 = numExtracted;		// the estimate number of key points cannot be greater than the number of available key points
    // compute a precise estimate keypoint number (num2)
    int num2 = 0;
    if (num1 > 0) {
        unsigned char tempbuf[MAX_DSC_LENGTH];
        BitOutputStream out2(tempbuf, MAX_DSC_LENGTH);			// attach temporary output buffer
        CsscCoordinateCoding cc(params);
        cc.generateHistogramMap(*this, num1);
        cc.toBinary(out2);			// write bitstream
        writeStream(out2, params.numRelevantPoints > 0, num1); // write features
        size_t numbits = out2.produced();						// ignore header bits
        double points = (availableBits * num1) / numbits;
        num2 = (int) (points + 0.25);		// round with a conservative approach (better less than more)
    }
    else
        throw CdvsException("No space left to encode local features");

    if (num2 > numExtracted)
        num2 = numExtracted;		// the estimate number of key points cannot be greater than the number of available key points

    if(params.selectMaxPoints > 0 && num2 > params.selectMaxPoints)
        return params.selectMaxPoints;

    return num2;
}