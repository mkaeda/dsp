#ifndef CONVOLVER_H
#define CONVOLVER_H

#include "wavfile.h"
#include <cstring>

class Convolver
{
public:
    Convolver();
    ~Convolver();
    void convolve(WavFile &x, WavFile &h, WavFile &y);
};

#endif
