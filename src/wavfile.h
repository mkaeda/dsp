#ifndef WAVFILE_H
#define WAVFILE_H

#include <fstream>
#include <string>
#include <iostream>
#include <cstring>

struct SubChunk1
{
    std::uint8_t chunk_id[4];
    std::uint32_t chunk_size;
    std::uint8_t format[4];
    std::uint8_t subchunk1_id[4];
    std::uint32_t subchunk1_size;
    std::uint16_t audio_format;
    std::uint16_t num_channels;
    std::uint32_t sample_rate;
    std::uint32_t byte_rate;
    std::uint16_t block_align;
    std::uint16_t bits_per_sample;
};

struct SubChunk2
{
    std::uint8_t subchunk2_id[4];
    std::uint32_t subchunk2_size;
};

class WavFile
{
public:
    WavFile(const char *filename);
    ~WavFile();
    void readHeader();
    void readData(short *buffer, int numSamples);
    void write(SubChunk1 sc1, SubChunk2 sc2, const short *data, int dataSize);

    SubChunk1 getHeaderSubChunk1() const;
    SubChunk2 getHeaderSubChunk2() const;
    int getNumSamples() const;

private:
    std::ifstream fileInputStream;
    std::string _filename;
    SubChunk1 subChunk1;
    SubChunk2 subChunk2;
    int numSamples;
};

#endif
