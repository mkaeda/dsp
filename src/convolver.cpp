#include "convolver.h"

#include <iostream>

Convolver::Convolver()
{
    // Constructor if needed
}

Convolver::~Convolver()
{
    // Destructor if needed
}

void convolveArrays(float x[], int N, float h[], int M, float y[], int P);
void convertDataToFloatArray(short *data, int dataSize, float *out, int outSize);
void convertFloatArrayToData(float *arr, int arrSize, short *out, int outSize);

void Convolver::convolve(WavFile &x, WavFile &h, WavFile &y)
{
    const int N = x.getNumSamples();
    const int M = h.getNumSamples();
    const int P = N + M - 1;

    short *xData = new short[N];
    short *hData = new short[M];
    short *yData = new short[P];

    float *xBuffer = new float[N];
    float *hBuffer = new float[M];
    float *yBuffer = new float[P];

    x.readData(xData, x.getHeaderSubChunk2().subchunk2_size);
    h.readData(hData, x.getHeaderSubChunk2().subchunk2_size);

    convertDataToFloatArray(xData, N, xBuffer, N);
    convertDataToFloatArray(hData, M, hBuffer, M);

    convolveArrays(xBuffer, N, hBuffer, M, yBuffer, P);

    // Update output header values.
    SubChunk1 sc1 = x.getHeaderSubChunk1();
    sc1.chunk_size = (P * sizeof(short)) + sizeof(SubChunk1) + sizeof(SubChunk2) - 8;
    sc1.subchunk1_size = 16;
    sc1.audio_format = 1;
    sc1.num_channels = 1;
    sc1.sample_rate = 44100;

    SubChunk2 sc2 = x.getHeaderSubChunk2();
    sc2.subchunk2_size = P * sizeof(short);

    convertFloatArrayToData(yBuffer, P, yData, P);

    y.write(sc1, sc2, yData, P * sizeof(short));

    delete[] xData;
    delete[] hData;
    delete[] yData;
    delete[] xBuffer;
    delete[] hBuffer;
    delete[] yBuffer;
}

void convertDataToFloatArray(short *data, int dataSize, float *out, int outSize)
{
    int bufSize = (outSize < dataSize) ? outSize : dataSize;
    for (int i = 0; i < bufSize; i++)
    {
        out[i] = static_cast<float>(data[i]) / 32768.0f;
    }
}

void convertFloatArrayToData(float *arr, int arrSize, short *out, int outSize)
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
            out[i] = static_cast<short>(arr[i] * 32768.0f);
        }
    }
}
void convolveArrays(float x[], int N, float h[], int M, float y[], int P)
{
    int n, m;

    /* Clear Output Buffer y[] */
    for (n = 0; n < P; n++)
    {
        y[n] = 0.0;
    }

    /* Outer Loop: process each input value x[n] in turn */
    for (n = 0; n < N; n++)
    {
        /* Inner loop: process x[n] with each sample of h[n] */
        for (m = 0; m < M; m++)
        {
            y[n + m] += x[n] * h[m];
        }
    }
}
