#include <iostream>
#include <vector>
#include <fstream>
#include "gnuplot-iostream.h"

struct WAVHeader {
    char chunkID[4];
    unsigned int chunkSize;
    char format[4];
    char subchunk1ID[4];
    unsigned int subchunk1Size;
    unsigned short audioFormat;
    unsigned short numChannels;
    unsigned int sampleRate;
    unsigned int byteRate;
    unsigned short blockAlign;
    unsigned short bitsPerSample;
    char subchunk2ID[4];
    unsigned int subchunk2Size;
};

// Hàm để đọc file WAV và trả về dữ liệu âm thanh
std::vector<double> readWAV(const std::string& filename, WAVHeader& header) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Couldn't open file: " << filename << std::endl;
        return {};
    }

    file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

    if (std::string(header.chunkID, 4) != "RIFF" ||
        std::string(header.format, 4) != "WAVE" ||
        std::string(header.subchunk1ID, 4) != "fmt " ||
        header.audioFormat != 1) {
        std::cerr << "Error: Invalid WAV file format." << std::endl;
        return {};
    }

    if (header.bitsPerSample != 16) {
        std::cerr << "Error: Only 16-bit audio is supported." << std::endl;
        return {};
    }

    std::vector<double> audioData;
    uint32_t numSamples = header.subchunk2Size / (header.bitsPerSample / 8);
    audioData.reserve(numSamples);

    int16_t sample;
    while (file.read(reinterpret_cast<char*>(&sample), sizeof(int16_t))) {
        audioData.push_back(static_cast<double>(sample) / std::pow(2, header.bitsPerSample - 1));
    }

    return audioData;
}

int main() {
    std::string filename = "D:\\PTIT\\lap_trinh_am_thanh\\ConsoleApplication1\\example1.wav";
    WAVHeader header;
    std::vector<double> audioData = readWAV(filename, header);

    Gnuplot gp;

    // Thiết lập đồ thị
    gp << "set title 'Waveform'\n";
    gp << "set xlabel 'Time (seconds)'\n";
    gp << "set ylabel 'Amplitude'\n";
    gp << "plot '-' with lines\n";

    double timeStep = 1.0 / header.sampleRate;
    double currentTime = 0.0;

    for (size_t i = 0; i < audioData.size(); ++i) {
        gp << currentTime << " " << audioData[i] << std::endl;
        currentTime += timeStep;
    }

    gp << "e\n";

    return 0;
}
