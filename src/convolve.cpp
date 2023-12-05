#include <iostream>
#include "wavfile.h"
#include "convolver.h"

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " sampleTone impulseTone\n";
        return -1;
    }

    WavFile sample(argv[1]);
    WavFile impulse(argv[2]);
    WavFile output(argv[3]);

    sample.readHeader();
    impulse.readHeader();

    Convolver convolver = Convolver();
    convolver.convolve(sample, impulse, output);

    return 0;
}