#pragma once

#include <vector>
#include "Parameters.h"
#include "Feature.h"
#include "FeatureListDB.h"
#include "PointPairs.h"
#include "BitOutputStream.h"
#include "BitInputStream.h"
#include <eigen3/Eigen/Dense>
#define HALF_PI (1.570796326794897)

class FeatureList
{
public:
	FeatureList(int targetLength = 0);
    void copyFeatureList(FeatureList& featureList);
    ~FeatureList();

    void setDescLength(size_t descLength);
	void setResolution(int imgWidth, int imgHeight, int originalMaxRes);
    void sortFeatures(bool (*compare)(const Feature&, const Feature&));
    static bool sortSpatialIndex(const Feature &, const Feature &);
    static bool sortRelevance(const Feature &, const Feature &);

	size_t NumExtracted() const;
    size_t size();
    void setRelevantPoints(int num);

    void append(const Feature& feature);
	void selectFirst(int n);

    void writeSiftStream(BitOutputStream &out) const;
    void readSiftStream(BitInputStream &in);
    void writeSift(const char * filename) const;
    void readSift(const char * filename);
    void writeStream(BitOutputStream &writer, bool writeRelevance, int nFeatures);
    void readStream(BitInputStream &reader, bool readRelevance);

	void compressHM(int targetLength);
    int  computeMaxPoints(const Parameters & params, int targetBits);

    int matchSift(PointPairs &pairs, FeatureList &otherFeatureList, float ratioThreshold);
	int matchOneWay(PointPairs &pairs, FeatureList &otherList, float ratioThreshold);
	int matchTwoWay(PointPairs &pairs, FeatureList &otherList, float ratioThreshold);

    int imageHeight;						///< the (possibly scaled) image height.
    int imageWidth;							///< the (possibly scaled) image width.
    int originalMaxResolution;				///< the maximum side resolution of the original image.
    vector<Feature> features;			    ///< the vector of SIFT features that describe the image.

private:
	friend class FeatureListDB;
    unsigned int qdescr_size;				///< The number of quantized elements in the key-point features (qdescr, H-mode)
    size_t numExtracted;					///< number of features already extracted from the source image
};