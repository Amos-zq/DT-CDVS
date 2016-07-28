/*
 * This software module was originally developed by:
 *
 *   Telecom Italia
 *
 * in the course of development of ISO/IEC <number> (Compact Descriptors for Visual 
 * Search) standard for reference purposes and its performance may not have been 
 * optimized. This software module includes implementation of one or more tools as 
 * specified by the ISO/IEC <number> standard.
 *
 * ISO/IEC gives you a royalty-free, worldwide, non-exclusive, copyright license to copy, 
 * distribute, and make derivative works of this software module or modifications thereof 
 * for use in implementations of the ISO/IEC <number> standard in products that satisfy 
 * conformance criteria (if any).
 *
 * Those intending to use this software module in products are advised that its use may 
 * infringe existing patents. ISO/IEC have no liability for use of this software module 
 * or modifications thereof.
 *
 * Copyright is not released for products that do not conform to audiovisual and image-
 * coding related ITU Recommendations and/or ISO/IEC International Standards.
 *
 ****** Section to be removed when the standard is published **************************
 *
 * Assurance that the originally developed software module can be used
 *  (1) in the ISO/IEC <number> standard once this standard has been adopted; and
 *  (2) to develop the ISO/IEC <number> standard:
 *
 * Telecom Italia grant ISO/IEC all rights necessary to include the originally 
 * developed software module or modifications thereof in the ISO/IEC <number> standard 
 * and to permit ISO/IEC to offer You a royalty-free, worldwide, non-exclusive, copyright 
 * license to copy, distribute, and make derivative works for use in implementations of 
 * the ISO/IEC <number> standard in products that satisfy conformance criteria (if any), 
 * and to the extent that such originally developed software module or portions of it 
 * are included in the ISO/IEC <number> standard.
 *
 * To the extent that Telecom Italia own patent rights that would be required 
 * to make, use, or sell the originally developed software module or portions thereof 
 * included in the ISO/IEC <number> standard in a conforming product, the <original 
 * developers> will assure the ISO/IEC that they are willing to negotiate licenses under 
 * reasonable and non-discriminatory terms and conditions with applicants throughout 
 * the world.
 *
 * ISO/IEC gives You a free license to this software module or modifications thereof 
 * for the sole purpose of developing the ISO/IEC <number> standard.
 *
 ****** End of section to be removed when the standard is published *******************
 *
 * Telecom Italia retain full rights to modify and use the code for their own 
 * purposes, assign or donate the code to a third party and to inhibit third parties 
 * from using the code for products that do not conform to MPEG-related 
 * ITU Recommendations and/or ISO/IEC International Standards.
 *
 * This copyright notice must be included in all copies or derivative works.
 * Copyright (c) ISO/IEC 2011.
 *
 */

#include "BitInputStream.h"
#include <string.h>
#include <assert.h>


BitInputStream::BitInputStream()
{
	isOpen				= false;		// intialize all member vars
	m_iBitBuffer			= 0;
	m_uNumOfBitsInBuffer	= 0;
	pBuf					= NULL;
	pStartBuf				= NULL;
	pEndBuf					= NULL;
}

BitInputStream::BitInputStream(const unsigned char * buf, size_t size)
{
	open(buf, size);
}

BitInputStream::~BitInputStream()
{
	if (isOpen)
		close();
}

void BitInputStream::open(const unsigned char * buf, size_t size)
{
	isOpen	= true;		// now is opened!
	m_iBitBuffer = 0;
	m_uNumOfBitsInBuffer = 0;
	pBuf = pStartBuf = buf;
	pEndBuf = buf + size;
}

void BitInputStream::close()
{
	m_iBitBuffer			= 0;
	m_uNumOfBitsInBuffer	= 0;
	pBuf					= NULL;
	pEndBuf					= NULL;
	isOpen = false;		// now is closed!
}

bool BitInputStream::eof () const
{
	assert(isOpen);
	if(m_uNumOfBitsInBuffer == 0)
		return (pBuf >= pEndBuf);
	else
		return false;
}

size_t BitInputStream::consumed() const
{
	return (8*(pBuf - pStartBuf) - m_uNumOfBitsInBuffer);
}

size_t BitInputStream::available() const
{
	return (m_uNumOfBitsInBuffer + 8*(pEndBuf - pBuf));
}

void BitInputStream::reset()
{
	assert(isOpen);
	m_iBitBuffer = 0;
	m_uNumOfBitsInBuffer = 0;
	pBuf = pStartBuf;
}

unsigned int BitInputStream::align()
{
	assert(isOpen);
	unsigned int numBits = m_uNumOfBitsInBuffer & 0x7;
	if (numBits > 0)
	{
		m_iBitBuffer <<= numBits;
		m_uNumOfBitsInBuffer -= numBits;
	}
	return numBits;
}

void BitInputStream::skip(unsigned int nbits)
{
	if (nbits > 0)
	{			
		if (nbits < m_uNumOfBitsInBuffer)		// the simple case...
		{
			m_iBitBuffer <<= nbits;
			m_uNumOfBitsInBuffer -= nbits;
		}
		else									// the difficult case...
		{
			jumpTo(consumed() + nbits);			// jump to the new absolute position
		}
	}
}


unsigned int BitInputStream::read()				// read one bit!
{
	assert(isOpen);
	assert(available() >= 1);			// cannot require more bits than those available!

	if(!m_uNumOfBitsInBuffer)
		fetch();						// if the bit buffer is empty, reload it
	int ret=(m_iBitBuffer >> 31);
	m_iBitBuffer <<= 1;
	m_uNumOfBitsInBuffer--;
	return ret;
}


unsigned int BitInputStream::read(unsigned int numBits)		// read n bits!
{
	assert(isOpen);
	assert(numBits <= available());		// cannot require more bits than those available!
	
	if ( !m_uNumOfBitsInBuffer)			// if the bit buffer is empty, reload it
		fetch();

	if(numBits <= m_uNumOfBitsInBuffer)
	{
		unsigned int tmp = m_iBitBuffer>>(32-numBits);
		m_iBitBuffer <<= numBits;
		m_uNumOfBitsInBuffer -= numBits;
		return tmp;
	}
	else
	{

		unsigned int tmp = m_iBitBuffer>>(32-m_uNumOfBitsInBuffer);
		numBits -= m_uNumOfBitsInBuffer;			// decrease numbits of what I have already read
		
		fetch();									// read next 32 bits
		
		tmp =(tmp<<numBits) | (m_iBitBuffer >> (32-numBits));
		m_iBitBuffer <<= numBits;
		m_uNumOfBitsInBuffer -= numBits;
		return tmp;
	}
}


void BitInputStream::read(unsigned char * destination, unsigned int nbits)		// block read!
{
	assert(isOpen);						// must be open
	assert(nbits <= available());		// cannot require more bits than those available!
	assert((nbits % 8) == 0);			// check if byte-aligned

	if (nbits >= 8)
	{
		memcpy(destination, getPointer(), nbits/8);
		skip(nbits);
	}
}


// fetch() - read bits form the bitstream in Big Endian format.
// this can only be called when the bitBuffer is empty.
// the method fetches 32 bits from the bitstream into the bit buffer. 
void BitInputStream::fetch()
{
	size_t bytes = pEndBuf - pBuf;		// currently available bytes
	if (bytes >= 4)
	{
		m_iBitBuffer = (pBuf[0] << 24) | (pBuf[1] << 16) | (pBuf[2] << 8) | pBuf[3]; 
		pBuf += 4;
		m_uNumOfBitsInBuffer = 32;
	}
	else	// in case there are less than 4 bytes to read, read what is available.
	{	
		m_iBitBuffer = 0;
		m_uNumOfBitsInBuffer = 0;

		while(pEndBuf > pBuf)
		{
			m_iBitBuffer = m_iBitBuffer << 8 | *(pBuf++);
			m_uNumOfBitsInBuffer += 8;
		}

		m_iBitBuffer <<= 8*(4 - bytes);
	}
}


const unsigned char * BitInputStream::getPointer() const
{
	assert((m_uNumOfBitsInBuffer % 8) == 0);				// test if byte-aligned 
	return (pBuf - (m_uNumOfBitsInBuffer / 8));				// return the current pointer
}

// jump to the absolute position (in bits) indicated in the param
void BitInputStream::jumpTo(size_t position)
{
	assert(isOpen);
	assert(position <= 8*getSize());

	unsigned int bitpos = position % 8 ;
	pBuf = pStartBuf + (position / 8);				// jump to the new position
	
	fetch();										// read from bitstream
	
	skip(bitpos);									// skip the remaining bits
}

size_t BitInputStream::getSize() const
{
	return pEndBuf - pStartBuf;
}
