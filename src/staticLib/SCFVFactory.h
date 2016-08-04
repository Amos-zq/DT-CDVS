#pragma once

#include "gaussian_mixture.h"
#include "FeatureList.h"
#include "SCFVSignature.h"

class SCFVFactory
{
private:
    static const int UNCOMPRESSED_SIFT_DESC_LENGTH = 128;
    static const int SCFV_MAX_NUMBER_FEATURES = 500;
    static const float GMM_W[];
    static const float GMM_MU[];
    static const float GMM_SIGMA[];
    static const float PCA_SIFT_POWER_LAW;
    static const float SCFV_POWER_LAW;
    static const float PCA_SIFT_MU[];
    static const float PCA_SIFT_EIGEN[];

    bool m_bHasVar;
    float m_fSCFVThreshold;
    gaussian_mixture mygmm;

    unsigned int generateSFV(const FeatureList& featureList, float* pFisherVector, bool vWordVisited[]);
public:
    SCFVFactory(const Parameters & params);		// constructor
    void generateSCFV(const FeatureList& featureList, SCFVSignature & signature);
};

