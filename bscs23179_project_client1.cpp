#include <iostream>
#include <winsock2.h>
#include <ws2tcpip.h>
#include <string>
#include <fstream>
#include <opencv2/opencv.hpp> // OpenCV header added

using namespace std;
#pragma comment(lib, "ws2_32.lib")

#define bufferCapacity 1024
#define sPortNo 9000

string getCurrentDateTime()
{
    time_t now = time(nullptr);
    tm localTime;
    localtime_s(&localTime, &now);

    char buffer[150];
    strftime(buffer, sizeof(buffer), "Date Recorded: %B %d, %Y - Time: %I:%M:%S %p", &localTime);

    return string(buffer);
}

void savePatientResult(string patientName, string inputDNA, string analysisResult)
{
    string fileName = "SampleResult_" + patientName + ".txt";
    ofstream outFile(fileName, ios::app);

    if (outFile.is_open())
    {
        outFile << "========== BioInformatics Result ==========" << endl;
        outFile << "Input DNA : " << inputDNA << endl;
        outFile << "DNA Analysis  : " << analysisResult << endl;
        outFile << getCurrentDateTime() << endl;
        outFile << "==========================================" << endl;
        outFile.close();
        cout << "\nResult saved in " << fileName << " successfully!\n";
    }
    else
    {
        cerr << "Error saving file!" << endl;
    }
}

// SIMPLIFIED VISUALIZATION FUNCTION
void visualizeSimpleDNA(const string& healthyDNA, const string& mutatedDNA, const string& patientName)
{
    int letterSpacing = 30; // Horizontal spacing between letters
    int margin = 60;        // Margin
    int imgWidth = max(healthyDNA.length(), mutatedDNA.length()) * letterSpacing + 100;
    int imgHeight = 120;

    cv::Mat vis(imgHeight, imgWidth, CV_8UC3, cv::Scalar(255, 255, 255)); // white background

    // Font settings
    int fontFace = cv::FONT_HERSHEY_SIMPLEX;
    double fontScale = 0.7;
    int thickness = 2;

    // Draw title
    cv::putText(vis, "DNA Mutation Visualization:", cv::Point(10, margin - 20),
        fontFace, fontScale, cv::Scalar(0, 0, 0), thickness);

    // Draw DNA letters in a single line
    for (int i = 0; i < healthyDNA.length(); i++)
    {
        cv::Scalar color;
        char displayChar;

        if (i < mutatedDNA.length() && healthyDNA[i] == mutatedDNA[i]) {
            // Match - green
            color = cv::Scalar(0, 200, 0);
            displayChar = mutatedDNA[i];
        }
        else {
            // Mismatch or gap - red
            color = cv::Scalar(0, 0, 255);
            displayChar = (i < mutatedDNA.length()) ? mutatedDNA[i] : '-';
        }

        int x = margin + i * letterSpacing;

        // Draw the character
        string letter(1, displayChar);
        cv::putText(vis, letter, cv::Point(x, margin + 20),
            fontFace, fontScale, color, thickness);
    }

    // Save and Show
    string outImg = "DNAVisualization_" + patientName + ".png";
    cv::imwrite(outImg, vis);
    cout << "DNA visualization saved " << outImg << endl;

    cv::imshow("DNA Visualization", vis);
    cv::waitKey(0);
    cv::destroyAllWindows();
}

int main()
{
    WSADATA wsa;
    SOCKET st;
    string inputG, patientName;
    struct sockaddr_in server;
    char bufferSpace[bufferCapacity] = { 0 };

    WSAStartup(MAKEWORD(2, 2), &wsa);
    st = socket(AF_INET, SOCK_STREAM, 0);

    server.sin_family = AF_INET;
    server.sin_port = htons(sPortNo);
    inet_pton(AF_INET, "127.0.0.1", &server.sin_addr);

    if (connect(st, (struct sockaddr*)&server, sizeof(server)) == SOCKET_ERROR)
    {
        cerr << "Connection couldn't establish." << endl;
        return 1;
    }

    cout << "Enter Sample Id: ";
    cin >> patientName;

    cout << "Enter DNA sequence: ";
    cin >> inputG;

    send(st, inputG.c_str(), inputG.length(), 0);
    int resultSize = recv(st, bufferSpace, bufferCapacity, 0);

    if (resultSize > 0)
    {
        bufferSpace[resultSize] = '\0';

        cout << "\n========== BioInformatics Result ==========" << endl;
        cout << "Input DNA : " << inputG << endl;
        cout << "DNA Analysis  : " << bufferSpace << endl;
        cout << "Date & Time : " << getCurrentDateTime() << endl;
        cout << "==========================================" << endl;

        savePatientResult(patientName, inputG, bufferSpace);

        // Call the simplified visualization
        visualizeSimpleDNA(bufferSpace, inputG, patientName);
    }
    else
    {
        cout << "No response from GeneLens System" << endl;
    }

    closesocket(st);
    WSACleanup();
    return 0;
}
