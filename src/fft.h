#ifndef CONVOLVER_H
#define CONVOLVER_H

#include "wavfile.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIZE 8
#define PI 3.141592653589793
#define TWO_PI (2.0 * PI)
#define SWAP(a, b) \
    tempr = (a);   \
    (a) = (b);     \
    (b) = tempr

class FftConvolver
{
public:
    FftConvolver();
    ~FftConvolver();
    void convolve(WavFile &x, WavFile &h, WavFile &y);
};

#endif
