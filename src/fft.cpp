#include "fft.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

FftConvolver::FftConvolver()
{
    // Constructor if needed
}

FftConvolver::~FftConvolver()
{
    // Destructor if needed
}

unsigned long getPaddedSize(int xSize, int hSize);
void fillComplex(short *data, int dataSize, std::vector<double> &complex);
void four1(double data[], int nn, int isign);

void FftConvolver::convolve(WavFile &x, WavFile &h, WavFile &y)
{
    short *xData = new short[x.getNumSamples()];
    short *hData = new short[h.getNumSamples()];

    x.readData(xData, x.getHeaderSubChunk2().subchunk2_size);
    h.readData(hData, h.getHeaderSubChunk2().subchunk2_size);

    const unsigned long N = getPaddedSize(x.getNumSamples(), h.getNumSamples());

    std::vector<double> xComplex(N * 2, 0.0);
    std::vector<double> hComplex(N * 2, 0.0);

    fillComplex(xData, x.getNumSamples(), xComplex);
    fillComplex(hData, h.getNumSamples(), hComplex);

    // Perform FFT on both complex arrays
    four1(&xComplex[0] - 1, N, 1);
    four1(&hComplex[0] - 1, N, 1);

    // Complex multiplication: result = x * h
    double real_x, imag_x, real_h, imag_h;
    for (unsigned long i = 0; i < N * 2; i += 2)
    {
        real_x = xComplex[i];
        imag_x = xComplex[i + 1];
        real_h = hComplex[i];
        imag_h = hComplex[i + 1];

        xComplex[i] = real_x * real_h - imag_x * imag_h;
        xComplex[i + 1] = real_x * imag_h + imag_x * real_h;
    }

    // Find the maximum absolute value in xComplex
    double maxAbsValue = 0.0;
    for (unsigned long i = 0; i < N * 2; ++i)
    {
        double absValue = std::abs(xComplex[i]);
        if (absValue > maxAbsValue)
        {
            maxAbsValue = absValue;
        }
    }

    // Normalize xComplex
    for (unsigned long i = 0; i < N * 2; ++i)
    {
        xComplex[i] /= maxAbsValue;
    }

    four1(&xComplex[0] - 1, N, -1);

    // After scaling
    for (unsigned long i = 0; i < N * 2; i += 2)
    {
        xComplex[i] = fmax(SHRT_MIN, fmin(SHRT_MAX, xComplex[i]));
        xComplex[i + 1] = fmax(SHRT_MIN, fmin(SHRT_MAX, xComplex[i + 1]));
    }

    // Before writing to the output file
    std::cout << "Min value in xComplex: " << *std::min_element(xComplex.begin(), xComplex.end()) << std::endl;
    std::cout << "Max value in xComplex: " << *std::max_element(xComplex.begin(), xComplex.end()) << std::endl;

    // Extract the real part of the result
    std::vector<short> result(N);
    for (unsigned long i = 0; i < N; ++i)
    {
        result[i] = static_cast<short>(round(xComplex[i * 2]));
    }

    // Before writing to the output file
    std::cout << "Min value in result: " << *std::min_element(result.begin(), result.end()) << std::endl;
    std::cout << "Max value in result: " << *std::max_element(result.begin(), result.end()) << std::endl;

    // Update output header values.
    SubChunk1 sc1 = x.getHeaderSubChunk1();
    sc1.chunk_size = (N * sizeof(short)) + sizeof(SubChunk1) + sizeof(SubChunk2) - 8;
    sc1.subchunk1_size = 16;
    sc1.audio_format = 1;
    sc1.num_channels = 1;
    sc1.sample_rate = 44100;

    SubChunk2 sc2 = x.getHeaderSubChunk2();
    sc2.subchunk2_size = (N * sizeof(short));

    y.write(sc1, sc2, &result[0], (N * sizeof(short)));

    delete[] xData;
    delete[] hData;
}

unsigned long getPaddedSize(int xSize, int hSize)
{
    int maxValue = std::max(xSize, hSize);
    unsigned long ret = 1;
    while (ret < (unsigned long)maxValue)
    {
        ret <<= 1;
    }
    return ret;
}

void fillComplex(short *data, int dataSize, std::vector<double> &complex)
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
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2)
    {
        if (j > i)
        {
            SWAP(data[j], data[i]);
            SWAP(data[j + 1], data[i + 1]);
        }
        m = nn;
        while (m >= 2 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while (n > mmax)
    {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2)
        {
            for (i = m; i <= n; i += istep)
            {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j + 1];
                tempi = wr * data[j + 1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j + 1] = data[i + 1] - tempi;
                data[i] += tempr;
                data[i + 1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}
