#include <iostream>
#include <vector>
#include <fstream>
#include <cmath> // For pow function
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
    const char* filename = "fntest.wav"; // Đường dẫn tới file WAV của bạn
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
    // In thông tin từ header
    cout << "______________READCURRENTFILE________________" << endl;
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
    std::cout << "Press ENTER to select MENU" << endl;
    cout << "---------------------------------------------" << endl;

    gp << "plot '-' with lines" << endl;
    gp.send(points);
    cin.get();
}
void readFile2() {
    const char* filename = "fntest2.wav"; // Đường dẫn tới file WAV của bạn
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
    // In thông tin từ header
    cout << "______________READCURRENTFILE________________" << endl;
    std::cout << "Number of channels: " << header.numChannels << std::endl;
    std::cout << "Sample rate: " << header.sampleRate << " Hz" << std::endl;
    std::cout << "Bits per sample: " << header.bitsPerSample << " bits" << std::endl;
    std::cout << "Duration: " << (double)header.subchunk2Size / (header.sampleRate * header.numChannels * header.bitsPerSample / 8) << " seconds" << std::endl;
    cout << "_____________________________________________" << endl;
    // Đọc dữ liệu âm thanh từ subchunk thứ hai của file WAV
    std::vector<int16_t> samples(header.subchunk2Size / sizeof(int16_t));
    file.read(reinterpret_cast<char*>(samples.data()), header.subchunk2Size);

    for (int i = 0; i < samples.size(); i++) {
        points2.push_back(make_pair(i, samples[i]));
    }
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

        gp << "plot '-' with lines" << endl;
        gp.send(tempPoints);
        cin.get();
    }

}

void sum(vector<pair<int, double>> vector1, vector<pair<int, double>> vector2) {
    formatVector(vector1, vector2);
    for (int i = 0; i < currentPointsAfterFormat.size(); i++) {
        points3.push_back(make_pair(currentPointsAfterFormat[i].first, currentPointsAfterFormat[i].second + pointsToMathAfterFormat[i].second));
    }
    gp << "plot '-' with lines" << endl;
    gp.send(points3);
    cin.get();
}

void minus(vector<pair<int, double>> vector1, vector<pair<int, double>> vector2) {
    formatVector(vector1, vector2);
    for (int i = 0; i < currentPointsAfterFormat.size(); i++) {
        points3.push_back(make_pair(currentPointsAfterFormat[i].first, currentPointsAfterFormat[i].second - pointsToMathAfterFormat[i].second));
    }
    gp << "plot '-' with lines" << endl;
    gp.send(points3);
    cin.get();
}

void times(vector<pair<int, double>> vector1, vector<pair<int, double>> vector2) {
    formatVector(vector1, vector2);
    for (int i = 0; i < currentPointsAfterFormat.size(); i++) {
        points3.push_back(make_pair(currentPointsAfterFormat[i].first, currentPointsAfterFormat[i].second * pointsToMathAfterFormat[i].second));
    }
    gp << "plot '-' with lines" << endl;
    gp.send(points3);
    cin.get();
}

int main() {
    readFile();
    readFile2();
    while (true) {
        std::cout << "_________________________MINI_PROJECT_____________________\n";
        std::cout << "              1. Original wav\n";
        std::cout << "              2. Time reversal\n";
        std::cout << "              3. Timeshift\n";
        std::cout << "              4. Down sampling\n";
        std::cout << "              5. Up sampling\n";
        std::cout << "              6. Sum\n";
        std::cout << "              7. Minus\n";
        std::cout << "              8. Times\n";
        std::cout << "              9. Quit\n";
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
            cout << "Enter time shift index" << endl;
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
            minus(points, points2);
            break;
        case 8:
            times(points, points2);
            break;
        case 9:
            break;
        default:
            break;
        }
        
    }
    return 0;
}

