#include <iostream>
#include <winsock2.h>
#include <ws2tcpip.h>
#include <string>
using namespace std;
#pragma comment(lib, "ws2_32.lib")


#define bufferCapacity 1024
#define sPortNo 9000

int main()
{
    cout << " ===============================\n";
    cout << "   Welcome to GeneLens !! \n";
    cout << " ===============================\n";

    WSADATA wsa;
    SOCKET st;
    string inputG;
    struct sockaddr_in server;
    char bufferSpace[bufferCapacity] = { 0 };

    WSAStartup(MAKEWORD(2, 2), &wsa);
    st = socket(AF_INET, SOCK_STREAM, 0);

    server.sin_family = AF_INET;
    server.sin_port = htons(sPortNo);
    inet_pton(AF_INET, "127.0.0.1", &server.sin_addr);

    if (connect(st, (struct sockaddr*)&server, sizeof(server)) == SOCKET_ERROR)
    {
        cerr << "connection couldnt establish." << endl;
        return 1;
    }

    cout << "Enter DNA sequence : ";
    cin >> inputG;

    send(st, inputG.c_str(), inputG.length(), 0);
    int resultSize = recv(st, bufferSpace, bufferCapacity, 0);

    if (resultSize > 0)
    {
        bufferSpace[resultSize] = '\0';

        cout << "DNA Analysis " << endl;
        cout << bufferSpace << endl;
    }
    else {
        cout << "No response from GeneLens System" << endl;
    }

    closesocket(st);
    WSACleanup();
    return 0;
}
