#include <iostream>
#include <winsock2.h>
#include <ws2tcpip.h>
#include <string>
#include <fstream>

using namespace std;
#pragma comment(lib, "ws2_32.lib")

#define bufferCapacity 1024
#define sPortNo 9001

string getCurrentDateTime()
{
    time_t now = time(nullptr);
    tm localTime;
    localtime_s(&localTime, &now);

    char buffer[150];
    strftime(buffer, sizeof(buffer), "Date Recorded: %B %d, %Y - Time: %I:%M:%S %p", &localTime);

    return string(buffer);
}

void saveResultFile(string sampleId, string inputDNA, string analysisResult)
{
    string fileName = "SampleResult_" + sampleId + ".txt";
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
        cout << "Input sample DNA : " << inputG << endl;
        cout << "DNA Analysis  : " << bufferSpace << endl;
        cout << "Date and Time : " << getCurrentDateTime() << endl;
        cout << "==========================================" << endl;

        saveResultFile(patientName, inputG, bufferSpace);
    }
    else
    {
        cout << "No response from GeneLens System" << endl;
    }

    closesocket(st);
    WSACleanup();
    return 0;
}
