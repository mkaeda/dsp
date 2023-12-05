#include "wavfile.h"

WavFile::WavFile(const char *filename)
{
    _filename = std::string(filename);
    fileInputStream.open(_filename, std::ios::binary);
}

WavFile::~WavFile()
{
    fileInputStream.close();
}

void WavFile::readHeader()
{
    fileInputStream.read(reinterpret_cast<char *>(&subChunk1), sizeof(SubChunk1));

    if (subChunk1.subchunk1_size != 16)
    {
        // eliminate null Bytes
        int remainder = subChunk1.subchunk1_size - 16;
        char randomVar[remainder];
        fileInputStream.read(randomVar, remainder);
    }

    fileInputStream.read(reinterpret_cast<char *>(&subChunk2), sizeof(SubChunk2));

    numSamples = subChunk2.subchunk2_size / (subChunk1.bits_per_sample / 8);
}

void WavFile::readData(short *buffer, int bufSize)
{
    std::ifstream file;
    file.open(_filename, std::ios::binary);
    // Read the audio data from the file
    file.read(reinterpret_cast<char *>(buffer), bufSize);
    file.close();
}

void WavFile::write(SubChunk1 sc1, SubChunk2 sc2, const short *data, int dataSize)
{
    if (static_cast<uint32_t>(dataSize) != sc2.subchunk2_size)
    {
        std::cerr << "Data size does not match file header.\n";
        return;
    }

    std::ofstream file;
    file.open(_filename, std::ios::binary);

    file.write(reinterpret_cast<const char *>(&sc1), sizeof(SubChunk1));

    if (sc1.subchunk1_size != 16)
    {
        std::cout << sc1.subchunk1_size << "\n";
        // write null bytes to it ends up as 16
        int remainder = sc1.subchunk1_size - 16;
        char randomVar[remainder];
        file.write(randomVar, remainder);
    }

    file.write(reinterpret_cast<const char *>(&sc2), sizeof(SubChunk2));
    file.write(reinterpret_cast<const char *>(data), dataSize);
    file.flush();
    file.close();
}

SubChunk1 WavFile::getHeaderSubChunk1() const
{
    return subChunk1;
}

SubChunk2 WavFile::getHeaderSubChunk2() const
{
    return subChunk2;
}

int WavFile::getNumSamples() const
{
    return numSamples;
}
