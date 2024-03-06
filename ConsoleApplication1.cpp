#include <iostream>
#include <vector>
#include <fstream>
#include <cmath> // For pow function
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

// Hàm để đảo ngược thứ tự của dữ liệu âm thanh
std::vector<double> reverseTime(const std::vector<double>& audioData) {
    std::vector<double> reversedData(audioData.rbegin(), audioData.rend());
    return reversedData;
}

std::vector<double> addDelay(const std::vector<double>& audioData, int delayMilliseconds, const WAVHeader& header) {
    std::vector<double> delayedData;
    delayedData.reserve(audioData.size());

    for (size_t i = 0; i < audioData.size(); ++i) {
        delayedData.push_back(audioData[i]);
    }

    // Thêm khoảng trễ vào đầu dữ liệu âm thanh
    int numSamplesToDelay = static_cast<int>((delayMilliseconds * 0.001) * header.sampleRate);
    for (int i = 0; i < numSamplesToDelay; ++i) {
        delayedData.insert(delayedData.begin(), 0); // Thêm giá trị 0 vào đầu
    }

    return delayedData;
}

// Hàm làm sớm thời gian trong mili giây
std::vector<double> shortenTime(const std::vector<double>& audioData, int timeAdvanceMilliseconds, const WAVHeader& header) {
    std::vector<double> shortenedData;
    shortenedData.reserve(audioData.size());

    // Tính số mẫu cần loại bỏ để làm sớm thời gian
    int numSamplesToRemove = static_cast<int>((timeAdvanceMilliseconds * 0.001) * header.sampleRate);

    // Copy dữ liệu âm thanh từ vị trí numSamplesToRemove đến cuối
    shortenedData.insert(shortenedData.end(), audioData.begin() + numSamplesToRemove, audioData.end());

    return shortenedData;
}


int main() {
    int number;
    std::string filename = "D:\\PTIT\\lap_trinh_am_thanh\\ConsoleApplication1\\example1.wav";
    WAVHeader header;
    std::vector<double> audioData = readWAV(filename, header);

    Gnuplot gp;

    std::cin >> number;

    if (number == 1) {
        // Thiết lập đồ thị
        gp << "set title 'Waveform'\n";
        gp << "set xlabel 'Time (seconds)'\n";
        gp << "set ylabel 'Amplitude'\n";
        gp << "plot '-' with lines\n";

        // Đảo ngược thứ tự của dữ liệu âm thanh
        std::vector<double> reversedData = reverseTime(audioData);

        double timeStep = 1.0 / header.sampleRate;
        double currentTime = 0.0; // Bắt đầu từ thời điểm đầu tiên

        for (size_t i = 0; i < reversedData.size(); ++i) {
            gp << currentTime << " " << reversedData[i] << std::endl;
            currentTime += timeStep; // Tăng thời gian
        }

        gp << "e\n";
    }
    else if (number == 2) {
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
    } else if (number == 3) {
        // Thêm trễ vào dữ liệu âm thanh
        std::vector<double> delayedData = addDelay(audioData, 1000, header); // Trễ 1000 miligiây, tương đương 1 giây

        // Thiết lập đồ thị
        gp << "set title 'Waveform'\n";
        gp << "set xlabel 'Time (seconds)'\n";
        gp << "set ylabel 'Amplitude'\n";
        gp << "plot '-' with lines\n";

        double timeStep = 1.0 / header.sampleRate;
        double currentTime = 0.0;

        for (size_t i = 0; i < delayedData.size(); ++i) {
            gp << currentTime << " " << delayedData[i] << std::endl;
            currentTime += timeStep;
        }

        gp << "e\n";
    } else if(number == 4) {
        std::vector<double> shortenedData = shortenTime(audioData, 1000, header); // Làm sớm 1000 miligiây, tương đương 1 giây

        // Thiết lập đồ thị
        gp << "set title 'Waveform'\n";
        gp << "set xlabel 'Time (seconds)'\n";
        gp << "set ylabel 'Amplitude'\n";
        gp << "plot '-' with lines\n";

        double timeStep = 1.0 / header.sampleRate;
        double currentTime = 0.0;

        for (size_t i = 0; i < shortenedData.size(); ++i) {
            gp << currentTime << " " << shortenedData[i] << std::endl;
            currentTime += timeStep;
        }

        gp << "e\n";
    }

    return 0;
}

