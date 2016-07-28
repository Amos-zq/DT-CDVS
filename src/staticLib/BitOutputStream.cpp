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

#include "BitOutputStream.h"
#include <assert.h>
#include <string.h>


BitOutputStream::BitOutputStream()
{
	isOpen = false;
	bitBuffer = 0; 
	availableBits = 32;
	startBuffer = endBuffer = writePtr = NULL;
}

BitOutputStream::BitOutputStream(unsigned char * buf, size_t size)
{
	open(buf, size);
}

BitOutputStream::~BitOutputStream()
{
	if (isOpen)
		close();
}

size_t BitOutputStream::produced() const
{
	return (writePtr - startBuffer)*8 + 32 - availableBits;
}

size_t BitOutputStream::available() const
{
	return (endBuffer - writePtr)*8 - 32 + availableBits;
}

void BitOutputStream::write(unsigned int data, unsigned int numBits)
{
	assert(numBits <= 32);			// must be contained in a 32 bit register
	assert(numBits > 0);			// must be at least 1
	assert((numBits == 32) ? true : ((data & ((unsigned int) 0xFFFFFFFF << numBits)) == 0));	// accept only data consistent with the bit number

	if(availableBits == 0)
		store();					// store 4 bytes into the write buffer

	if(numBits <= availableBits)			// the simple case
	{
		availableBits -= numBits;
		bitBuffer = (bitBuffer << numBits) | data;
		if (availableBits == 0)
			store();
	}
	else									// the (numBits > available) case
	{
		numBits -= availableBits;
		bitBuffer = (bitBuffer << availableBits) | (data >> numBits) ;
		store();							// in this case store() is called with wrong "availableBits" - do not put "asserts" in there.
		availableBits = 32 - numBits;
		bitBuffer = data;
	}
}

void BitOutputStream::write(unsigned int bit)		// write one bit!
{
	assert(isOpen);
	assert((bit==0)||(bit==1));

	if(availableBits == 0)
		store();					// store 4 bytes into the write buffer

	availableBits--;
	bitBuffer = (bitBuffer << 1) | bit;
}

void BitOutputStream::write(unsigned char * source, unsigned int nbits)
{
	assert(isOpen);						// must be open
	assert(nbits <= available());		// cannot require more bits than those available!
	assert((nbits % 8) == 0);			// check if byte-aligned

	if (nbits >= 8)
	{
		memcpy(getPointer(), source , nbits/8);
		skip(nbits);
	}
}

void BitOutputStream::reset()
{
	assert(isOpen);
	availableBits = 32;
	bitBuffer = 0;
	writePtr = startBuffer;			// but keep the byte buffer
}

void BitOutputStream::open(unsigned char * buf, size_t size)
{
	isOpen = true;
	startBuffer = buf;
	endBuffer = buf + size;
	reset();
}

size_t BitOutputStream::close()
{
	flush(0);
	size_t producedBits = produced();
	bitBuffer = 0; 
	availableBits = 32;
	writePtr = startBuffer = NULL;
	isOpen = false;
	return producedBits;
}

bool BitOutputStream::eof () const
{
	assert(isOpen);

	return ((writePtr + (32 - availableBits) / 8) >= endBuffer);
}

void BitOutputStream::align(unsigned int fill)
{
	unsigned int nBits = availableBits & 7;
	if (nBits) 
	{
		if (fill)
			write((1<<(nBits-1))-1, nBits);
		else
			write((unsigned int) 0, nBits);
	}
}

// this can only be called when the bitBuffer is full.
// the method stores 32 bits into the bitstream buffer. 
void BitOutputStream::store()
{
	// ignore availableBits here - the bit buffer is full in any case!

	size_t bytes = endBuffer - writePtr;			// currently available bytes
	assert(bytes >= 4);								// cannot write more than expected!
	if (bytes >= 4)
	{
		// at least 4 bytes of the write buffer are free

		*(writePtr++) = (unsigned char) (bitBuffer >> 24);
		*(writePtr++) = (unsigned char) (bitBuffer >> 16);
		*(writePtr++) = (unsigned char) (bitBuffer >> 8);
		*(writePtr++) = (unsigned char) bitBuffer ;
	}
	
	bitBuffer = 0;
	availableBits = 32;	
}

// this can be called at any time, including at the end (e.g. by close).
void BitOutputStream::flush(unsigned int fill)
{
	assert(isOpen);
	
	align(fill);										// align to the next byte boundary

	unsigned char * tmpPtr = writePtr;					// do NOT change the write pointer!
	unsigned int writtenBits = 32 - availableBits;

	for (unsigned int i = 0; i < writtenBits; i += 8)
		*(tmpPtr++) = (unsigned char) (bitBuffer >> ((writtenBits - 8 - i)));
}

unsigned char * BitOutputStream::getPointer() const
{
	assert((availableBits % 8) == 0);				// test if byte-aligned 
	return writePtr + ((32 - availableBits) / 8);	// return the current pointer
}

// jump to the absolute position (in bits) indicated in the param
void BitOutputStream::jumpTo(size_t position)
{
	assert(isOpen);
	assert(position <= 8*getSize());
	assert((position % 8) == 0);					// must be byte aligned

	flush(0);										// flush any remaining bit
	writePtr = startBuffer + (position / 8);		// jump to the new position
	bitBuffer = 0;									// reset the write register
	availableBits = 32;		
}

size_t BitOutputStream::getSize() const
{
	return endBuffer - startBuffer;
}



void BitOutputStream::skip(unsigned int nbits)
{
	if (nbits > 0)
	{	
		size_t newAddress = produced() + nbits;		// new bit address
		assert((newAddress % 8) == 0);
		jumpTo(newAddress);									// jump to the new absolute position
	}
}
