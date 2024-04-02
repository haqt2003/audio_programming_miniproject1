#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <utility>
#include <numeric>
#include "gnuplot-iostream.h"

using namespace std;

Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");
std::vector<pair<int, double>>points;
std::vector<pair<int, double>>points2;
std::vector<pair<int, double>>points3;
vector<pair<int, double>>tempPoints;
vector<pair<int, double>>tempPointsOfVec1;
vector<pair<int, double>>tempPointsOfVec2;
vector<pair<int, double>>currentPointsAfterFormat;
vector<pair<int, double>>pointsToMathAfterFormat;
double samplingFrequency;
unsigned int sRate;
unsigned short nChannels;
unsigned short bPerSample;

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

void readFile() {
    const char* filename = "original.wav"; // Đường dẫn tới file WAV của bạn
    std::ifstream file(filename, std::ios::binary);

    if (!file) {
        std::cerr << "Error: Failed to open file " << filename << std::endl;
    }
    // Đọc header của file WAV
    WAVHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

    // Kiểm tra định dạng WAV
    if (std::string(header.chunkID, 4) != "RIFF" || std::string(header.format, 4) != "WAVE") {
        std::cerr << "Error: Unsupported WAV format" << std::endl;
    }
    samplingFrequency = header.sampleRate;
    sRate = header.sampleRate;
    nChannels = header.numChannels;
    bPerSample = header.bitsPerSample;
    // In thông tin từ header
    cout << "______________" << filename << "________________" << endl;
    std::cout << "Number of channels: " << header.numChannels << std::endl;
    std::cout << "Sample rate: " << header.sampleRate << " Hz" << std::endl;
    std::cout << "Bits per sample: " << header.bitsPerSample << " bits" << std::endl;
    std::cout << "Duration: " << (double)header.subchunk2Size / (header.sampleRate * header.numChannels * header.bitsPerSample / 8) << " seconds" << std::endl;
    cout << "_____________________________________________" << endl;
    // Đọc dữ liệu âm thanh từ subchunk thứ hai của file WAV
    std::vector<int16_t> samples(header.subchunk2Size / sizeof(int16_t));
    file.read(reinterpret_cast<char*>(samples.data()), header.subchunk2Size);

    for (int i = 0; i < samples.size(); i++) {
        points.push_back(make_pair(i, samples[i]));
    }

    //test
    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Drawing..." << endl;
    std::cout << "Press ENTER to select MENU" << endl;

    gp << "plot '-' with lines" << endl;
    gp.send(points);
    cin.get();
}

void readFile2() {
    const char* filename = "original2.wav"; // Đường dẫn tới file WAV của bạn
    std::ifstream file(filename, std::ios::binary);

    if (!file) {
        std::cerr << "Error: Failed to open file " << filename << std::endl;
    }
    // Đọc header của file WAV
    WAVHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

    // Kiểm tra định dạng WAV
    if (std::string(header.chunkID, 4) != "RIFF" || std::string(header.format, 4) != "WAVE") {
        std::cerr << "Error: Unsupported WAV format" << std::endl;
    }
    
    // Đọc dữ liệu âm thanh từ subchunk thứ hai của file WAV
    std::vector<int16_t> samples(header.subchunk2Size / sizeof(int16_t));
    file.read(reinterpret_cast<char*>(samples.data()), header.subchunk2Size);

    for (int i = 0; i < samples.size(); i++) {
        points2.push_back(make_pair(i, samples[i]));
    }
}

void writeFile(const std::string& filename, const std::vector<int16_t>& samples, unsigned int sampleRate, unsigned short numChannels, unsigned short bitsPerSample) {
    // Mở file để ghi dữ liệu
    std::ofstream file(filename, std::ios::binary);

    if (!file) {
        std::cerr << "Error: Failed to open file for writing: " << filename << std::endl;
        return;
    }

    // Tính toán kích thước của dữ liệu âm thanh
    unsigned int dataSize = samples.size() * sizeof(int16_t);

    // Tạo header cho file WAV
    WAVHeader header;
    memcpy(header.chunkID, "RIFF", 4);
    header.chunkSize = dataSize + 36;
    memcpy(header.format, "WAVE", 4);
    memcpy(header.subchunk1ID, "fmt ", 4);
    header.subchunk1Size = 16;
    header.audioFormat = 1;
    header.numChannels = numChannels;
    header.sampleRate = sampleRate;
    header.bitsPerSample = bitsPerSample;
    header.byteRate = sampleRate * numChannels * bitsPerSample / 8;
    header.blockAlign = numChannels * bitsPerSample / 8;
    memcpy(header.subchunk2ID, "data", 4);
    header.subchunk2Size = dataSize;

    // Ghi header vào file
    file.write(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

    // Ghi dữ liệu âm thanh vào file
    file.write(reinterpret_cast<const char*>(samples.data()), dataSize);

    // Đóng file sau khi ghi xong
    file.close();

    std::cout << "Successfully wrote WAV file: " << filename << std::endl;
}

int findRoot(const vector<pair<int, double>>& nameV) {
    int tmp = 0;
    for (int i = 0; i < nameV.size(); i++) {
        if (nameV[i].first == 0) {
            tmp = i;
            break;
        }
    }
    return tmp;
}

void formatVector(vector<pair<int, double>> vector1, vector<pair<int, double>> vector2) {

    tempPointsOfVec1.clear();
    tempPointsOfVec2.clear();
    tempPoints.clear();

    int posOfVct1 = findRoot(vector1); //default_wav_file
    int posOfVct2 = findRoot(vector2);
    int coutOfNumberVct1 = vector1.size() - posOfVct1 - 1;
    int coutOfNumberVct2 = vector2.size() - posOfVct2 - 1;

    if (posOfVct1 < posOfVct2) {
        for (int i = vector2[0].first; i < vector1[0].first; i++) {
            tempPointsOfVec1.push_back(make_pair(i, 0));
        }
        for (int i = 0; i < vector1.size(); i++) {
            tempPointsOfVec1.push_back(make_pair(vector1[i].first, vector1[i].second));
        }
        vector1.clear();
        vector1 = tempPointsOfVec1;

    }
    else if (posOfVct1 > posOfVct2) {
        for (int i = vector1[0].first; i < vector2[0].first; i++) {
            tempPointsOfVec2.push_back(make_pair(i, 0));
        }
        for (int i = 0; i < vector2.size(); i++) {
            tempPointsOfVec2.push_back(make_pair(vector1[i].first, vector2[i].second));
        }
        vector2.clear();
        vector2 = tempPointsOfVec2;
    }

    if (coutOfNumberVct1 < coutOfNumberVct2) {
        tempPointsOfVec1.clear();
        for (int i = 0; i < vector1.size(); i++) {
            tempPointsOfVec1.push_back(make_pair(vector1[i].first, vector1[i].second));
        }

        for (int i = tempPointsOfVec1[tempPointsOfVec1.size() - 1].first + 1; i <= vector2[vector2.size() - 1].first; i++) {
            tempPointsOfVec1.push_back(make_pair(i, 0));
        }
        vector1.clear();
        vector1 = tempPointsOfVec1;
    }
    else if (coutOfNumberVct1 > coutOfNumberVct2) {
        tempPointsOfVec2.clear();

        for (int i = 0; i < vector2.size(); i++) {
            tempPointsOfVec2.push_back(make_pair(vector2[i].first, vector2[i].second));
        }

        for (int i = tempPointsOfVec2[tempPointsOfVec2.size() - 1].first + 1; i <= vector1[vector1.size() - 1].first; i++) {
            tempPointsOfVec2.push_back(make_pair(i, 0));
        }
        vector2.clear();
        vector2 = tempPointsOfVec2;
    }
    currentPointsAfterFormat.clear();
    tempPoints = vector2;
    currentPointsAfterFormat = tempPoints;
    tempPoints.clear();

    pointsToMathAfterFormat.clear();
    tempPoints = vector1;
    pointsToMathAfterFormat = tempPoints;
    tempPoints.clear();
}

void reverseFile() {
    tempPoints.clear();
    for (int i = 0; i < points.size(); i++) {
        tempPoints.push_back(make_pair(-points[i].first, points[i].second));
        sort(tempPoints.begin(), tempPoints.begin() + i + 1);
    }
    vector<int16_t> writeVec;
    for (int i = 0; i < tempPoints.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(tempPoints[i].second));
    }
    writeFile("Reverse.wav", writeVec, sRate, nChannels, bPerSample);
    points.clear();
    points = tempPoints;

    gp << "plot '-' with lines" << endl;
    gp.send(tempPoints);
    cin.get();
}

void timeShift(int n) {
    tempPoints.clear();
    int pos;
    int cnt;

    if (n <= 0) {
        n = abs(n);
        pos = findRoot(points);
        cnt = points.size() - 1 - pos;
        if(n <= cnt) {
            for (int i = 0; i < points.size(); i++) {
                tempPoints.push_back(make_pair(points[i].first - n, points[i].second));
            }
            points.clear();
            points = tempPoints;
            gp << "plot '-' with lines" << endl;
            gp.send(tempPoints);
            cin.get();
        }
        else if(n > cnt) {
            for (int i = 0; i < points.size(); i++) {
                tempPoints.push_back(make_pair(points[i].first, points[i].second));
            }
            int addNum = n - cnt;
            while (addNum > 0) {
                tempPoints.push_back(make_pair(0, 0));
                addNum--;
            }
            int vl = 0;
            for (int i = tempPoints.size() - 1; i >= 0; i--) {
                tempPoints[i].first = vl;
                vl--;
            }
            points.clear();
            points = tempPoints;
            gp << "plot '-' with lines" << endl;
            gp.send(tempPoints);
            cin.get();
        }
    }
    else {
        n = abs(n);
        pos = findRoot(points);
        cnt = pos - 0;
        if (n <= cnt) {
            for (int i = 0; i < points.size(); i++) {
                tempPoints.push_back(make_pair(points[i].first + n, points[i].second));
            }
            points.clear();
            points = tempPoints;
            gp << "plot '-' with lines" << endl;
            gp.send(tempPoints);
            cin.get();
        }
        else if (n > cnt) {
            int addNum = n - cnt;
            while (addNum > 0) {
                tempPoints.push_back(make_pair(0, 0));
                addNum--;
            }
            for (int i = 0; i < points.size(); i++) {
                tempPoints.push_back(make_pair(points[i].first, points[i].second));
            }
            int vl = 0;
            for (int i = 0; i < tempPoints.size(); i++) {
                tempPoints[i].first = vl;
                vl++;
            }
            vector<int16_t> writeVec;
            for (int i = 0; i < tempPoints.size(); i++) {
                writeVec.push_back(static_cast<int16_t>(tempPoints[i].second));
            }
            writeFile("Timeshift.wav", writeVec, sRate, nChannels, bPerSample);
            points.clear();
            points.clear();
            points = tempPoints;
            gp << "plot '-' with lines" << endl;
            gp.send(tempPoints);
            cin.get();
        }
    }
}

void downSample(int n) {
    tempPoints.clear();
    int pos = findRoot(points);

    for (int i = pos; i >= 0; i-=n) {
        tempPoints.push_back(make_pair(points[i].first / n, points[i].second));
    }
    reverse(tempPoints.begin(), tempPoints.end());
    for (int i = pos + n; i < points.size(); i += n) {
        tempPoints.push_back(make_pair(points[i].first / n, points[i].second));
    }
    vector<int16_t> writeVec;
    for (int i = 0; i < tempPoints.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(tempPoints[i].second));
    }
    writeFile("Downsample.wav", writeVec, sRate, nChannels, bPerSample);
    points.clear();
    points.clear();
    points = tempPoints;

    gp << "plot '-' with lines" << endl;
    gp.send(tempPoints);
    cin.get();
}

void upSample(int n) {
    tempPoints.clear();

    int pos = findRoot(points);
    int cnt = points.size() - pos - 1;

    if (pos > 0) {
        for (int i = 0; i < pos; i++) {
            tempPoints.push_back(make_pair(points[i].first * n, points[i].second));
            for (int j = 1; j < n; j++) {
                tempPoints.push_back(make_pair(points[i].first * n + j, 0));
            }
        }
        tempPoints.push_back(make_pair(0, points[pos].second));
        for (int i = 1; i <= cnt; i++) {
            for (int j = 1; j < n; j++) {
                tempPoints.push_back(make_pair(tempPoints[tempPoints.size() - 1].first + 1, 0));
            }
            tempPoints.push_back(make_pair(i * n, points[pos + i].second));
        }
        vector<int16_t> writeVec;
        for (int i = 0; i < tempPoints.size(); i++) {
            writeVec.push_back(static_cast<int16_t>(tempPoints[i].second));
        }
        writeFile("Upsample.wav", writeVec, sRate, nChannels, bPerSample);
        points.clear();
        gp << "plot '-' with lines" << endl;
        gp.send(tempPoints);
        cin.get();

    }
    else {
        tempPoints.push_back(make_pair(0, points[pos].second));
        for (int i = 1; i <= cnt; i++) {
            for (int j = 1; j < n; j++) {
                tempPoints.push_back(make_pair(tempPoints[tempPoints.size() - 1].first + 1, 0));
            }
            tempPoints.push_back(make_pair(i * n, points[pos + i].second));
        }
        vector<int16_t> writeVec;
        for (int i = 0; i < tempPoints.size(); i++) {
            writeVec.push_back(static_cast<int16_t>(tempPoints[i].second));
        }
        writeFile("Upsample.wav", writeVec, sRate, nChannels, bPerSample);
        points.clear();
        gp << "plot '-' with lines" << endl;
        gp.send(tempPoints);
        cin.get();
    }

}

void sum(vector<pair<int, double>> vector1, vector<pair<int, double>> vector2) {
    points3.clear();
    formatVector(vector1, vector2);
    for (int i = 0; i < currentPointsAfterFormat.size(); i++) {
        points3.push_back(make_pair(currentPointsAfterFormat[i].first, currentPointsAfterFormat[i].second + pointsToMathAfterFormat[i].second));
    }
    vector<int16_t> writeVec;
    for (int i = 0; i < points3.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(points3[i].second));
    }
    writeFile("Sum.wav", writeVec, sRate, nChannels, bPerSample);
    points.clear();
    gp << "plot '-' with lines" << endl;
    gp.send(points3);
    cin.get();
}

void times(vector<pair<int, double>> vector1, vector<pair<int, double>> vector2) {
    points3.clear();
    formatVector(vector1, vector2);
    for (int i = 0; i < currentPointsAfterFormat.size(); i++) {
        points3.push_back(make_pair(currentPointsAfterFormat[i].first, currentPointsAfterFormat[i].second * pointsToMathAfterFormat[i].second));
    }
    vector<int16_t> writeVec;
    for (int i = 0; i < points3.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(points3[i].second));
    }
    writeFile("Times.wav", writeVec, sRate, nChannels, bPerSample);
    points.clear();
    gp << "plot '-' with lines" << endl;
    gp.send(points3);
    cin.get();
}

void amplifySignal(double amplificationFactor) {
    tempPoints.clear();
    tempPoints = points;
    for (auto& sample : tempPoints) {
        sample.second *= amplificationFactor;
    }
    vector<int16_t> writeVec;
    for (int i = 0; i < tempPoints.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(tempPoints[i].second));
    }
    writeFile("Amplify Signal.wav", writeVec, sRate, nChannels, bPerSample);
    points.clear();
    points = tempPoints;
    gp << "plot '-' with lines" << endl;
    gp.send(tempPoints);
    cin.get();
}

vector <double> convolution(const vector <double>& v1,const vector<double>& v2) {
    vector<double> tmp(v1.size() + v2.size() - 1, 0);

    for (int i = 0; i < v1.size(); i++) {
        for (int j = 0; j < v2.size(); j++) {
            tmp[i + j] += v1[i] * v2[j];
        }
    }
    
    return tmp;
}

void LPF(int n, double cutSample) {
    tempPoints.clear();

    vector <double> filter;
    for (int i = 0; i < n; i++) {
        double wc = 2 * 3.1416 * cutSample / samplingFrequency;
        double h = (wc / 3.1416) * (sin(wc * (i - (n - 1) / 2.0)) / (wc * (i - (n - 1) / 2.0)));
        filter.push_back(h);
    }

    // Nhân chập
    vector <double> op(points.size() + filter.size() - 1, 0);
    for (int i = 0; i < filter.size(); i++) {
        for (int j = 0; j < points.size(); j++) {
            op[i+j] += filter[i] * points[j].second;
        }
    }
    
    // Ghi file
    vector<int16_t> writeVec;
    for (int i = 0; i < op.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(op[i]));
    }
    writeFile("LPF.wav", writeVec, sRate, nChannels, bPerSample);
    points.clear();
    points = tempPoints;

    // Vẽ đồ thị
    gp << "plot '-' with lines" << endl;
    gp.send(op);
    cin.get();
}

void HPF(int n, double cutSample) {
    tempPoints.clear();
    vector <double> filter;
    for (int i = 0; i < n; i++) {
        double wc = 2 * 3.1416 * cutSample / samplingFrequency;
        double omg;
        if (i - (n - 1) / 2 == 0) {
            omg = 1;
        }
        else {
            omg = 0;
        }
        double h = omg - ((wc / 3.1416) * ((sin(wc * (i - (n - 1) / 2))) / (wc * (i - (n - 1) / 2))));
        filter.push_back(h);
    }

    vector <double> op(points.size() + filter.size() - 1, 0);
    for (int i = 0; i < filter.size(); i++) {
        for (int j = 0; j < points.size(); j++) {
            op[i + j] += filter[i] * points[j].second;
        }
    }
    vector<int16_t> writeVec;
    for (int i = 0; i < op.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(op[i]));
    }
    writeFile("HPF.wav", writeVec, sRate, nChannels, bPerSample);
    points.clear();
    points = tempPoints;

    gp << "plot '-' with lines" << endl;
    gp.send(op);
    cin.get();
}

void BPF(int n, double cutSample1, double cutSample2) {
    tempPoints.clear();
    vector <double> filter;
    for (int i = 0; i < n; i++) {
        double wc1 = 2 * 3.1416 * cutSample1 / samplingFrequency;
        double wc2 = 2 * 3.1416 * cutSample2 / samplingFrequency;
        double h  = (wc2 / 3.1416) * ((sin(wc2 * (i - (n - 1) / 2))) / (wc2 * (i - (n - 1) / 2))) - (wc1 / 3.1416) * ((sin(wc1 * (i - (n - 1) / 2))) / (wc1 * (i - (n - 1) / 2)));
        filter.push_back(h);
    }
    vector <double> op(points.size() + filter.size() - 1, 0);
    for (int i = 0; i < filter.size(); i++) {
        for (int j = 0; j < points.size(); j++) {
            op[i + j] += filter[i] * points[j].second;
        }
    }
    vector<int16_t> writeVec;
    for (int i = 0; i < op.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(op[i]));
    }
    writeFile("BPF.wav", writeVec, sRate, nChannels, bPerSample);

    points.clear();
    points = tempPoints;

    gp << "plot '-' with lines" << endl;
    gp.send(op);
    cin.get();
}

void BSF(int n, double cutSample1, double cutSample2) {
    tempPoints.clear();
    vector <double> filter;
    for (int i = 0; i < n; i++) {
        double wc1 = 2 * 3.1416 * cutSample1 / samplingFrequency;
        double wc2 = 2 * 3.1416 * cutSample2 / samplingFrequency;
        double omg;
        if (i - (n - 1) / 2 == 0) {
            omg = 1;
        }
        else {
            omg = 0;
        }
        double h = omg - ((wc2 / 3.1416) * ((sin(wc2 * (i - (n - 1) / 2))) / (wc2 * (i - (n - 1) / 2))) - (wc1 / 3.1416) * ((sin(wc1 * (i - (n - 1) / 2))) / (wc1 * (i - (n - 1) / 2))));
        filter.push_back(h);
    }
    vector <double> op(points.size() + filter.size() - 1, 0);
    for (int i = 0; i < filter.size(); i++) {
        for (int j = 0; j < points.size(); j++) {
            op[i + j] += filter[i] * points[j].second;
        }
    }
    vector<int16_t> writeVec;
    for (int i = 0; i < op.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(op[i]));
    }
    writeFile("BSF.wav", writeVec, sRate, nChannels, bPerSample);

    points.clear();
    points = tempPoints;

    gp << "plot '-' with lines" << endl;
    gp.send(op);
    cin.get();
}

void echoEffect(int delaySamples, double decayFactor) {
    // Tạo vector tạm thời để lưu âm thanh sau khi áp dụng hiệu ứng echo
    vector<pair<int, double>> tempAudio(points.size() + delaySamples);

    // Áp dụng hiệu ứng echo
    for (size_t i = 0; i < points.size(); ++i) {
        // Copy âm thanh gốc vào vector tạm thời
        tempAudio[i].first = points[i].first;
        tempAudio[i].second = points[i].second;

        // Áp dụng echo cho âm thanh
        if (i >= delaySamples) {
            tempAudio[i].second += decayFactor * points[i - delaySamples].second;
        }
    }

    points.clear();
    tempPoints.clear();
    tempPoints = tempAudio;
    points = tempPoints;
    vector<int16_t> writeVec;
    for (int i = 0; i < tempPoints.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(tempPoints[i].second));
    }
    writeFile("Echo.wav", writeVec, sRate, nChannels, bPerSample);
    gp << "plot '-' with lines" << endl;
    gp.send(tempPoints);
    cin.get();
}

void reverbEffect(int delaySamples, double decayFactor) {
    // Tạo vector tạm thời để lưu âm thanh sau khi áp dụng hiệu ứng reverb
    vector<pair<int, double>> tempAudio(points.size() + delaySamples);

    // Áp dụng hiệu ứng reverb
    for (size_t i = 0; i < points.size(); ++i) {
        // Copy âm thanh gốc vào vector tạm thời
        tempAudio[i].first = points[i].first;
        tempAudio[i].second = points[i].second;

        // Áp dụng reverb cho âm thanh
        if (i >= delaySamples) {
            // Tính toán giá trị của âm thanh sau khi reverb
            double reverbValue = 0.0;
            for (int j = 0; j < delaySamples; ++j) {
                // Áp dụng hệ số suy giảm cho âm thanh đã phát
                reverbValue += decayFactor * points[i - delaySamples + j].second;
            }
            // Thêm âm thanh sau khi reverb vào vector tạm thời
            tempAudio[i].second += reverbValue;
        }
    }

    points.clear();
    tempPoints.clear();
    tempPoints = tempAudio;
    points = tempPoints;
    vector<int16_t> writeVec;
    for (int i = 0; i < tempPoints.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(tempPoints[i].second));
    }
    writeFile("Reverb.wav", writeVec, sRate, nChannels, bPerSample);
    gp << "plot '-' with lines" << endl;
    gp.send(tempPoints);
    cin.get();
}

void fadeIn(int fadeInDuration) {
    tempPoints.clear();
    // Xác định số lượng mẫu âm thanh
    int numSamples = points.size();

    // Nếu khoảng thời gian fade in lớn hơn số lượng mẫu, chỉnh lại để tránh truy cập ngoài phạm vi
    fadeInDuration = min(fadeInDuration, numSamples);

    tempPoints = points;
    // Duyệt qua các mẫu âm thanh và tăng dần âm lượng từ giá trị ban đầu lên giá trị tối đa
    for (int i = 0; i < fadeInDuration; ++i) {
        double ratio = static_cast<double>(i) / fadeInDuration;
        tempPoints[i].second *= ratio;
    }
    vector<int16_t> writeVec;
    for (int i = 0; i < tempPoints.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(tempPoints[i].second));
    }
    writeFile("FadeIn.wav", writeVec, sRate, nChannels, bPerSample);
    points.clear();
    points = tempPoints;

    gp << "plot '-' with lines" << endl;
    gp.send(tempPoints);
    cin.get();
}

void fadeOut(int fadeOutDuration) {
    tempPoints.clear();

    int numSamples = points.size();
    fadeOutDuration = min(numSamples, fadeOutDuration);

    tempPoints = points;

    for (int i = numSamples - fadeOutDuration; i < tempPoints.size(); i++) {
        double ratio = static_cast<double>(numSamples - i) / fadeOutDuration;
        tempPoints[i].second *= ratio;
    }
    tempPoints[tempPoints.size() - 1].second = 0;
    vector<int16_t> writeVec;
    for (int i = 0; i < tempPoints.size(); i++) {
        writeVec.push_back(static_cast<int16_t>(tempPoints[i].second));
    }
    writeFile("FadeOut.wav", writeVec, sRate, nChannels, bPerSample);
    points.clear();
    points = tempPoints;

    gp << "plot '-' with lines" << endl;
    gp.send(tempPoints);
    cin.get();
}

int main() {
    readFile();
    readFile2();
    while (true) {
        std::cout << "_________________________MINI PROJECT 1_____________________\n";
        std::cout << "              1. Original wav\n";
        std::cout << "              2. Time reverse\n";
        std::cout << "              3. Timeshift\n";
        std::cout << "              4. Down sampling\n";
        std::cout << "              5. Up sampling\n";
        std::cout << "              6. Sum\n";
        std::cout << "              7. Times\n";
        std::cout << "              8. Amplify signal\n";
        std::cout << "              9. LPF\n";
        std::cout << "              10. HPF\n";
        std::cout << "              11. BSF\n";
        std::cout << "              12. BPF\n";
        std::cout << "              13. Echo\n";
        std::cout << "              14. Reverb\n";
        std::cout << "              15. Fade in\n";
        std::cout << "              16. Fade out\n";
        std::cout << "              17. Quit\n";
        std::cout << "__________________________________________________________\n";
        std::cout << "Select: ";

        int select;
        cin >> select;

        switch (select)
        {
        case 1:
            readFile();
            break;
        case 2:
            reverseFile();
            break;
        case 3:
            cout << "Enter time shift sample:" << endl;
            int n;
            cin >> n;
            timeShift(n);
            break;
        case 4:
            cout << "Enter value:" << endl;
            int m;
            cin >> m;
            downSample(m);
            break;
        case 5:
            cout << "Enter value:" << endl;
            int l;
            cin >> l;
            upSample(l);
            break;
        case 6:
            sum(points, points2);
            break;
        case 7:
            times(points, points2);
            break;
        case 8:
            double k;
            cin >> k;
            amplifySignal(k);
            break;
        case 9:
            cout << "Enter filterOrder and cutoffFrequency value:" << endl;
            int n1;
            double cutSample;
            cin >> n1 >> cutSample;
            LPF(n1, cutSample);
            break;
        case 10:
            cout << "Enter filterOrder and cutoffFrequency value:" << endl;
            int n2;
            double cutSample2;
            cin >> n2 >> cutSample2;
            HPF(n2, cutSample2);
            break;
        case 11:
            cout << "Enter filterOrder, cutoffFrequency1 and cutoffFrequency1 value:" << endl;
            int n3;
            double cutSample31;
            double cutSample32;
            cin >> n3>> cutSample31 >> cutSample32;
            BSF(n3, cutSample31, cutSample32);
            break;
        case 12:
            cout << "Enter filterOrder, cutoffFrequency1 and cutoffFrequency1 value:" << endl;
            int n4;
            double cutSample41;
            double cutSample42;
            cin >> n4 >> cutSample41 >> cutSample42;
            BPF(n4, cutSample41, cutSample42);
            break;
        case 13:
            cout << "Enter delay and decay value:" << endl;
            int dl;
            double dc;
            cin >> dl >> dc;
            echoEffect(dl, dc);
            break;
        case 14:
            cout << "Enter delay and decay value:" << endl;
            int dl1;
            double dc1;
            cin >> dl1 >> dc1;
            reverbEffect(dl1, dc1);
            break;
        case 15:
            cout << "Enter duration value:" << endl;
            int duration;
            cin >> duration;
            fadeIn(duration);
            break;
        case 16:
            cout << "Enter duration value:" << endl;
            int duration1;
            cin >> duration1;
            fadeOut(duration1);
            break;
        case 17:
            break;
        default:
            break;
        }
        
    }
    return 0;
}

