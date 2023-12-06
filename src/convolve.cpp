#include <iostream>
#include "wavfile.h"
#include "fft.h"

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cerr << "usage: convolve inputfile IRfile outputfile\n";
        return -1;
    }

    WavFile sample(argv[1]);
    WavFile impulse(argv[2]);
    WavFile output(argv[3]);

    sample.readHeader();
    impulse.readHeader();

    FftConvolver convolver = FftConvolver();
    convolver.convolve(sample, impulse, output);

    return 0;
}