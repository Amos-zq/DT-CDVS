#pragma once
#include "defines.h"
#include "BitInputStream.h"
#include "BitOutputStream.h"

#define AC_MAX_VALUE4 81

class FeatureList;
class FeatureListDB;

class Feature {
public:
	Feature(void);

	void toFile(FILE * file) const;
	void fromFile(FILE * file);
    void toStream(BitOutputStream& out) const;
	void fromStream(BitInputStream& out);

	static const unsigned int descrLength = 128;		///< the size of a feature (key point)
	float x;				  ///< the X coordinate of the SIFT point
	float y;				  ///< the Y coordinate of the SIFT point
	float scale;			  ///< the scale of the SIFT point
	float orientation;		  ///< the orientation of the SIFT point
	float peak;				  ///< the peak of the SIFT point
	float curvRatio;		  ///< the ratio of the curvatures
#ifdef ADD_ATTRIBUTE
	float log_hessian_eign1;  ///< the eign value 1 of the LoG hessian
	float log_hessian_eign2;  ///< the eign value 2 of the LoG hessian
	float g_hessian_eign1;    ///< the eign value 2 of the Gaussian hessian
	float g_hessian_eign2;    ///< the eign value 2 of the Gaussian hessian
	float log_hessian_det;    ///< the determinant value of the LoG hessian
	float log_hessian_trace;  ///< the trace value of the LoG hessian
	float g_hessian_det;      ///< the determinant value of the Gaussian hessian
	float g_hessian_trace;    ///< the trace value of the Gaussian hessian
	float log_dx;
	float log_dy;
	float log_ds;
	float log_dxx;
	float log_dyy;
	float log_dss;
	float log_dxy;
	float log_dxs;
	float log_dys;
	float g_dx;
	float g_dy;
	float g_dxx;
	float g_dyy;
	float g_dxy;
#endif
	float descr[descrLength]; ///< the descriptor of the SIFT point
	float pdf;				  ///< probability of this point to be matched
	int spatialIndex;		  ///< indicates the order of transmission of this point
	unsigned short relevance; ///< relevance of the key-point, computed on the basis of his characteristics
	unsigned char ucdescr[32];			///<  the quantized descriptor used for storing and matching. 4 elements per byte.
	int qdescr[HM_SIFT_QUANT_SIZE];		///<  the quantized (ternary) descriptor values.
	int octave;				///< octave of this feature
	int iscale;				///< int scale
private:
	friend class FeatureList;
};


class FeatureDB {
private:
    friend class FeatureList;
    friend class FeatureListDB;
    static const unsigned short nQuantElements = 10;
    unsigned char quantElements[nQuantElements];      ///< quantized ternary elements of the SIFT descriptor.
public:
    unsigned short x;	///< x coordinate of the point.
    unsigned short y;	///< y coordinate of the point.
};
