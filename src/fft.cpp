#include "fft.h"

#include <iostream>
#include <algorithm>

FftConvolver::FftConvolver()
{
    // Constructor if needed
}

FftConvolver::~FftConvolver()
{
    // Destructor if needed
}

bool isPowerOfTwo(int n);
int roundToNextPowerOfTwo(int n);
void fillComplex(short *data, int dataSize, double *complex);
void convertComplexToData(double *complex, int arrSize, short *out, int outSize);
void four1(double data[], int nn, int isign);

void FftConvolver::convolve(WavFile &x, WavFile &h, WavFile &y)
{
    const int N = x.getNumSamples();
    const int M = h.getNumSamples();
    const int P = N + M - 1;
    const unsigned int pSize = roundToNextPowerOfTwo(std::max(N, M)) * 2;

    std::cout << "N=" << N << " M=" << M << " pSize=" << pSize << std::endl;

    short *xData = new short[N];
    short *hData = new short[M];
    x.readData(xData, N);
    h.readData(hData, M);

    double *xComplex = new double[pSize];
    double *hComplex = new double[pSize];

    fillComplex(xData, N, xComplex);
    fillComplex(hData, M, hComplex);

    // Perform FFT on both arrays
    four1(xComplex - 1, pSize, 1);
    // four1(hComplex-1, pSize, 1);

    // for (unsigned int i = 0; i < pSize; i += 2)
    // {
    //     // Perform element-wise multiplication on real and imaginary parts
    //     double real = xComplex[i] * hComplex[i];
    //     double imag = xComplex[i] * hComplex[i + 1];
    //     xComplex[i] = real / pSize;
    //     xComplex[i + 1] = imag / pSize;
    // }

    // double maxVal = *std::max_element(xComplex, xComplex + pSize);
    // double minVal = *std::min_element(xComplex, xComplex + pSize);
    // std::cout << "Max Value: " << maxVal << ", Min Value: " << minVal << std::endl;

    // // Inverse FFT on the multiplied array
    // four1(xComplex, pSize, -1);

    // // Update output header values.
    // SubChunk1 sc1 = x.getHeaderSubChunk1();
    // sc1.chunk_size = pSize + sizeof(SubChunk1) + sizeof(SubChunk2) - 8;
    // sc1.subchunk1_size = 16;
    // sc1.audio_format = 1;
    // sc1.num_channels = 1;
    // sc1.sample_rate = 44100;

    // SubChunk2 sc2 = x.getHeaderSubChunk2();
    // sc2.subchunk2_size = pSize;

    // y.write(sc1, sc2, (short *)xComplex, pSize);

    delete[] xData;
    delete[] hData;
    delete[] xComplex;
    delete[] hComplex;
}

bool isPowerOfTwo(int n)
{
    // Check if the integer is greater than 0 and has only one bit set
    return (n > 0) && ((n & (n - 1)) == 0);
}

int roundToNextPowerOfTwo(int n)
{
    int power = ceil(std::log2(n));
    return pow(2, power);
}

void fillComplex(short *data, int dataSize, double *complex)
{
    // Fill the complex arrays with the real data
    for (int i = 0; i < dataSize; i++)
    {
        complex[i * 2] = data[i];
        complex[i * 2 + 1] = 0.0;
    }
}

void four1(double data[], int nn, int isign)
{
    if (std::log2(nn) != floor(std::log2(nn)))
    {
        std::cout << "nn must be a power of 2\n";
        return;
    }
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2)
    {
        // if (j > i)
        // {
        //     SWAP(data[j], data[i]);
        //     SWAP(data[j + 1], data[i + 1]);
        // }
        // m = nn;
        // while (m >= 2 && j > m)
        // {
        //     j -= m;
        //     m >>= 1;
        // }
        // j += m;
    }
    std::cout << "i=" << i << "\n";

    // mmax = 2;
    // while (n > mmax)
    // {
    //     istep = mmax << 1;
    //     theta = isign * (6.28318530717959 / mmax);
    //     wtemp = sin(0.5 * theta);
    //     wpr = -2.0 * wtemp * wtemp;
    //     wpi = sin(theta);
    //     wr = 1.0;
    //     wi = 0.0;
    //     for (m = 1; m < mmax; m += 2)
    //     {
    //         for (i = m; i <= n; i += istep)
    //         {
    //             j = i + mmax;
    //             tempr = wr * data[j] - wi * data[j + 1];
    //             tempi = wr * data[j + 1] + wi * data[j];
    //             data[j] = data[i] - tempr;
    //             data[j + 1] = data[i + 1] - tempi;
    //             data[i] += tempr;
    //             data[i + 1] += tempi;
    //         }
    //         wr = (wtemp = wr) * wpr - wi * wpi + wr;
    //         wi = wi * wpr + wtemp * wpi + wi;
    //     }
    //     mmax = istep;
    // }
}
