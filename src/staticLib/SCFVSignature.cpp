#include "defines.h"
#include "SCFVSignature.h"
#include "LookupTable.h"
extern LookupTable lookupTable;
// ====================== class SCFVSignature ======================//

SCFVSignature::SCFVSignature():bHasVar(false) {
    clear();			// clear all data
}

SCFVSignature::SCFVSignature(bool myHasVar) {
    clear();			// clear all data
    hasVar(myHasVar);	// set the correct hasVar() value
}

void SCFVSignature::copySignature(SCFVSignature& signature) {
    bHasVar = signature.bHasVar;
    fNorm   = signature.fNorm;
    m_numVisited = signature.m_numVisited;
    for(int i = 0; i < numberCentroids; ++i) {
        m_vWordBlock[i] = signature.m_vWordBlock[i];
        m_vWordVarBlock[i] = signature.m_vWordVarBlock[i];
    }
}

size_t SCFVSignature::getVisited() const {
    return m_numVisited;
}

void SCFVSignature::setVisited() {
    int	count = 0;
    for (size_t k=0; k < numberCentroids; ++k) {
        if (m_vWordBlock[k])
            ++count;
    }
    m_numVisited = count;
}

float SCFVSignature::getNorm() const {
    return fNorm;
}

void SCFVSignature::setNorm() {
    setVisited();											// compute visited
    fNorm = sqrt((float) (getVisited() * PCASiftLength));	// compute norm
}

size_t SCFVSignature::size() const {
    if (hasVar())
        return 2 * numberCentroids * PCASiftLength;
    else
        return numberCentroids * PCASiftLength;
}

void SCFVSignature::clear() {
    m_numVisited = 0;
    fNorm = 0.0f;
    memset(m_vWordBlock, 0, sizeof(m_vWordBlock));				// reset all values to zero
    memset(m_vWordVarBlock, 0, sizeof(m_vWordVarBlock));		// reset all values to zero
}

bool SCFVSignature::hasVar()  const {
    return bHasVar;
}

void SCFVSignature::hasVar(bool value) {
    bHasVar = value;
}

void SCFVSignature::write(BitOutputStream & out) const {
    for (int k=0; k < numberCentroids; ++k)
        out.write(m_vWordBlock[k] > 0 ? 1 : 0);		// write visited bit

    for (int k=0; k < numberCentroids; ++k) {
        if (m_vWordBlock[k])
            out.write(m_vWordBlock[k], PCASiftLength);
    }

    if (hasVar()) {
        for (int k=0; k < numberCentroids; ++k) {
            if (m_vWordBlock[k])
                out.write(m_vWordVarBlock[k], PCASiftLength);
        }
    }
}

void SCFVSignature::read(BitInputStream & in) {
    bool visited[numberCentroids];			// temporary array
    memset(visited, 0, sizeof(visited));	// clear all values

    for (size_t k=0; k < numberCentroids; ++k)
        visited[k] = (in.read() == 1);		// read visited bits

    for (int k=0; k < numberCentroids; ++k) {
        if (visited[k])
            m_vWordBlock[k] = in.read(PCASiftLength);
    }

    if (hasVar()) {
        for (int k=0; k < numberCentroids; ++k) {
            if (visited[k])
                m_vWordVarBlock[k] = in.read(PCASiftLength);
        }
    }
    setNorm();		// compute the norm of this signature
}

void SCFVSignature::writeDB(FILE * file) const {
    fwrite(&bHasVar, sizeof(bHasVar), 1, file);
    fwrite(&fNorm, sizeof(fNorm), 1, file);
    fwrite(&m_numVisited, sizeof(m_numVisited), 1, file);
    fwrite(&m_vWordBlock, sizeof(m_vWordBlock), 1, file);
    fwrite(&m_vWordVarBlock, sizeof(m_vWordVarBlock), 1, file);
}

void SCFVSignature::readDB(FILE * file) {
    fread(&bHasVar, sizeof(bHasVar), 1, file);
    fread(&fNorm, sizeof(fNorm), 1, file);
    fread(&m_numVisited, sizeof(m_numVisited), 1, file);
    fread(&m_vWordBlock, sizeof(m_vWordBlock), 1, file);
    fread(&m_vWordVarBlock, sizeof(m_vWordVarBlock), 1, file);
}

void SCFVSignature::writeDB(BitOutputStream & writer) const {
	writer.write((unsigned char*)&bHasVar, sizeof(bHasVar) * 8);
	writer.write((unsigned char*)&fNorm, sizeof(fNorm) * 8);
	writer.write((unsigned char*)&m_numVisited, sizeof(m_numVisited) * 8);
	writer.write((unsigned char*)&m_vWordBlock, sizeof(m_vWordBlock) * 8);
	writer.write((unsigned char*)&m_vWordVarBlock, sizeof(m_vWordVarBlock) * 8);
}
void SCFVSignature::readDB(BitInputStream & reader) {
	reader.read((unsigned char*)&bHasVar, sizeof(bHasVar) * 8);
	reader.read((unsigned char*)&fNorm, sizeof(fNorm) * 8);
	reader.read((unsigned char*)&m_numVisited, sizeof(m_numVisited) * 8);
	reader.read((unsigned char*)&m_vWordBlock, sizeof(m_vWordBlock) * 8);
	reader.read((unsigned char*)&m_vWordVarBlock, sizeof(m_vWordVarBlock) * 8);
}
size_t SCFVSignature::sizeDB() const{
	int size = 0;
	size += sizeof(bHasVar);
	size += sizeof(fNorm);
	size += sizeof(m_numVisited);
	size += sizeof(m_vWordBlock);
	size += sizeof(m_vWordVarBlock);
	return size;
}


float SCFVSignature::compare(const SCFVSignature& referSignature) const{
    if(referSignature.getVisited() <= 5) return 2.5;

    bool bothHasVar = bHasVar && referSignature.hasVar();
    float totalCorrelation = 0.0, totalVarCorrelation = 0.0;
    unsigned int queryWord, referWord, xorWord, hanmingDist;

    for(int nCentroid = 0; nCentroid < numberCentroids; ++nCentroid) {
        queryWord = m_vWordBlock[nCentroid];
        referWord = referSignature.m_vWordBlock[nCentroid];
        if(queryWord && referWord) {
            xorWord   = queryWord ^ referWord;
            hanmingDist = lookupTable.fHanming[xorWord >> 16] + lookupTable.fHanming[xorWord & 65535];
            totalCorrelation += lookupTable.fCorrTable[hanmingDist];
            if(bothHasVar) {
                queryWord = m_vWordVarBlock[nCentroid];
                referWord = referSignature.m_vWordVarBlock[nCentroid];
                xorWord   = queryWord ^ referWord;
                hanmingDist = lookupTable.fHanming[xorWord >> 16] + lookupTable.fHanming[xorWord & 65535];
                totalVarCorrelation += lookupTable.fVarCorrTable[hanmingDist];
            }
        }
    }
    totalCorrelation /= (fNorm * referSignature.getNorm());
    totalVarCorrelation /= (fNorm * referSignature.getNorm());
    if(bothHasVar)  return (2 - (totalCorrelation + totalVarCorrelation));
    else            return (2 - 2* totalCorrelation);
}

float SCFVSignature::predictScore(size_t sumVisited, Parameters* params1, Parameters* params2, bool isMixed){
    if(isMixed){
#ifdef USE_PARA_REDUCTION
        int s = shift((params1->scfvThreshold + params2->scfvThreshold)/2, LookupTable::shiftParams[0], LookupTable::shiftParams[1]);
        return fourier(sumVisited + s, LookupTable::fourierParams[0] + params1->shiftOffsetMixed, LookupTable::fourierParams[1], LookupTable::fourierParams[2], LookupTable::fourierParams[3]);
#else
        return fourier(sumVisited, params1.matchParamsMixed[0], params2.matchParamsMixed[1], params1.matchParamsMixed[2], params2.matchParamsMixed[3]);
#endif
    } else {
#ifdef USE_PARA_REDUCTION
        float shiftX = params1->descLength == 512 ? 0 : params1->scfvThreshold;
        int s = shift(shiftX, LookupTable::shiftParams[0], LookupTable::shiftParams[1]);
        return fourier(sumVisited + s, LookupTable::fourierParams[0] + params1->shiftOffsetMixed, LookupTable::fourierParams[1], LookupTable::fourierParams[2], LookupTable::fourierParams[3]);
#else
        return fourier(sumVisited, params1.matchParams[0], params2.matchParams[1], params1.matchParams[2], params2.matchParams[3]);
#endif
    }
}

float SCFVSignature::fourier(float x, float a0, float a1, float b1, float w) {
    return a0 + a1 * cos(x * w) + b1 * sin(x * w);
}
int  SCFVSignature::shift(float mode, float a, float b) {
    if (mode <= 0)  return 0;
    else return -int(a * exp(b * mode) + 0.5);
}