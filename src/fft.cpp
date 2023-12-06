#include "fft.h"

#include <iostream>

FftConvolver::FftConvolver()
{
    // Constructor if needed
}

FftConvolver::~FftConvolver()
{
    // Destructor if needed
}

bool isPowerOfTwo(int n);
unsigned int roundUpToNextPowerOfTwo(unsigned int n);
void fillComplex(short *data, int dataSize, double *complex);
// void convertDataToDoubleArray(short *data, int dataSize, double *out, int outSize);
void convertDoubleArrayToData(double *arr, int arrSize, short *out, int outSize);
void createComplexSawtooth(double data[], int size);
void four1(double data[], int nn, int isign);

void FftConvolver::convolve(WavFile &x, WavFile &h, WavFile &y)
{
    const int N = isPowerOfTwo(x.getHeaderSubChunk2().subchunk2_size)
                      ? x.getHeaderSubChunk2().subchunk2_size
                      : roundUpToNextPowerOfTwo(x.getHeaderSubChunk2().subchunk2_size);
    const int M = isPowerOfTwo(h.getHeaderSubChunk2().subchunk2_size)
                      ? h.getHeaderSubChunk2().subchunk2_size
                      : roundUpToNextPowerOfTwo(h.getHeaderSubChunk2().subchunk2_size);
    const int P = isPowerOfTwo(N + M - 1)
                      ? N + M - 1
                      : roundUpToNextPowerOfTwo(N + M - 1);

    short *xData = new short[N];
    short *hData = new short[M];
    short *yData = new short[P * 2];

    std::fill(xData, xData + N, 0);
    std::fill(hData, hData + M, 0);

    x.readData(xData, N);
    h.readData(hData, M);

    double *xBuffer = new double[P * 2];
    double *hBuffer = new double[P * 2];

    fillComplex(xData, N, xBuffer);
    fillComplex(hData, M, hBuffer);
    // convertDataToDoubleArray(xData, N, xBuffer, N);
    // convertDataToDoubleArray(hData, M, hBuffer, M);

    // Perform FFT on both arrays
    four1(xBuffer, P, 1);
    four1(hBuffer, P, 1);

    // Element-wise multiplication in the frequency domain
    for (int i = 0; i < P * 2; i += 2)
    {
        double real = xBuffer[i] * hBuffer[i] - xBuffer[i + 1] * hBuffer[i + 1];
        double imag = xBuffer[i] * hBuffer[i + 1] + xBuffer[i + 1] * hBuffer[i];
        xBuffer[i] = real;
        xBuffer[i + 1] = imag;
    }

    // Inverse FFT on the multiplied array
    four1(xBuffer, P, -1);

    convertDoubleArrayToData(xBuffer, P * 2, yData, P * 2);

    // Update output header values.
    SubChunk1 sc1 = x.getHeaderSubChunk1();
    sc1.chunk_size = P * 2 + sizeof(SubChunk1) + sizeof(SubChunk2) - 8;
    sc1.subchunk1_size = 16;
    sc1.audio_format = 1;
    sc1.num_channels = 1;
    sc1.sample_rate = 44100;

    SubChunk2 sc2 = x.getHeaderSubChunk2();
    sc2.subchunk2_size = P * 2;

    y.write(sc1, sc2, yData, P * 2);

    delete[] xData;
    delete[] hData;
    delete[] yData;
    delete[] xBuffer;
    delete[] hBuffer;
}

bool isPowerOfTwo(int n)
{
    // Check if the integer is greater than 0 and has only one bit set
    return (n > 0) && ((n & (n - 1)) == 0);
}

unsigned int roundUpToNextPowerOfTwo(unsigned int n)
{
    if (n == 0)
    {
        return 1; // 2^0 = 1
    }

    n--; // Ensure that if n is already a power of 2, we don't need to add 1
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;

    return n + 1;
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

void convertDataToDoubleArray(short *data, int dataSize, double *out, int outSize)
{
    int bufSize = (outSize < dataSize) ? outSize : dataSize;
    for (int i = 0; i < bufSize; i++)
    {
        out[i] = (double)data[i];
    }
}

void convertDoubleArrayToData(double *arr, int arrSize, short *out, int outSize)
{
    int bufSize = (outSize < arrSize) ? outSize : arrSize;
    for (int i = 0; i < bufSize; i++)
    {
        if (arr[i] > 1)
        {
            out[i] = 32767;
        }
        else if (arr[i] < -1)
        {
            out[i] = -32678;
        }
        else
        {
            out[i] = static_cast<short>(arr[i] * 32768.0);
        }
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
