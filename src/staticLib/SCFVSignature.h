#pragma once
#include "Parameters.h"
#include "BitOutputStream.h"
#include "BitInputStream.h"

class SCFVSignature {
public:
    static const int PCASiftLength = 32;
#ifdef USE_NEW_GLOBAL
    static const int numberCentroids = 256;
#else
    static const int numberCentroids = 128; 
#endif		
    bool bHasVar;
    float fNorm;
    unsigned int m_numVisited;
    unsigned int m_vWordBlock[numberCentroids];
    unsigned int m_vWordVarBlock[numberCentroids];

    SCFVSignature();	// must be private to force using the constructor with parameters
    SCFVSignature(bool hasVar);	///< constructor which declares if this signature contains variance information
    void copySignature(SCFVSignature& signature);
    void clear();				///< clear all data
    size_t size() const;		///< get the size of the binary signature
    size_t getVisited() const;		///< get the number of visited words 
    void setVisited();						///< compute and store the correct number of visited words
    float getNorm() const;				///< get the norm of this signature
    void  setNorm();					///< compute and store the correct norm for this signature
    bool hasVar() const;				///< tell if this signature has variance information
    void hasVar(bool value);			///< set this signature as one containing variance information (if value is true) 
    void write(BitOutputStream & out) const;	///< write the binary signature into the given output stream
    void read(BitInputStream & in);				///< read the binary signature from the given input stream
    void writeDB(FILE * file) const;			///< write the signature to file
    void readDB(FILE * file);				///< read the signature from file
	void writeDB(BitOutputStream & writer) const;	///< write the binary signature into the given output stream
	void readDB(BitInputStream & reader);				///< read the binary signature from the given input stream
	size_t sizeDB() const;		///< get the size of the binary signature

    float compare(const SCFVSignature& referSignature) const;
    static float predictScore(size_t sumVisited, Parameters* params1, Parameters* params2, bool isMixed);
    static float fourier(float x, float a0, float a1, float b1, float w);
    static int   shift(float mode, float a, float b);
};