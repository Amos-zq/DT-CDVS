#pragma once
#include "BitInputStream.h"
#include "BitOutputStream.h"
#include "FeatureList.h"
#include "PointPairs.h"

class FeatureListDB
{
public:
    static const double pPos;
    static const double pNeg;
    static const double pRonPos[100];
    static const double pDonPos[100];
    static const double pRonNeg[100];
    static const double pDonNeg[100];
    static const double rdWeight[100 * 30];

    static const double dWeight[100];
    static const double rWeight[100];
    FeatureListDB(const char* filename = NULL);
    FeatureListDB(const FeatureListDB&);
    ~FeatureListDB(void);

    void copyFeatures(const FeatureList& a);
    void setFilename(const char *filename);
    int  nFeatures() const;

    int matchOneWay(const FeatureListDB &otherFeatureDBList, float ratioThreshold, PointPairs& pointPair);
    int matchTwoWay(const FeatureListDB &otherFeatureDBList, float ratioThreshold, PointPairs& pointPair);
    double expMatch(const FeatureListDB &otherFeatureDBList, PointPairs& pointPair, PointPairs& expPair, int nTest, float threshold);
    int findFeature(float x, float y) const;

	void write(FILE* file) const;
	void read(FILE* file);
	void write(BitOutputStream& writer) const;
	void read(BitInputStream& reader);
	int  size() const;

    std::string imagefile;				///< pathname of the image file.
    int numFeatures;					///< number of features of this image
    int nDescLength;					///< descriptor length in bytes.
    FeatureDB *features;				///< vector of stored features that describe this image.
};
