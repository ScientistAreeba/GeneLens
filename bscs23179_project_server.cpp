#include <iostream>
#include <winsock2.h>
#include <ws2tcpip.h>
#include <thread>
#include <vector>
using namespace std;

#pragma comment(lib, "ws2_32.lib")
#define bufferCapacity 1024

void forwardToWorker(const string& inputGene, SOCKET clientSocket, int preferredWorkerPort)
{
    int fallbackPort = 8000;

    // Try preferred port first
    SOCKET workerSocket = socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in workerAddr{};
    workerAddr.sin_family = AF_INET;
    workerAddr.sin_port = htons(preferredWorkerPort);
    inet_pton(AF_INET, "127.0.0.1", &workerAddr.sin_addr);

    if (connect(workerSocket, (sockaddr*)&workerAddr, sizeof(workerAddr)) == SOCKET_ERROR)
    {
        cerr << "Server: Worker on port " << preferredWorkerPort << " not available. Trying fallback port " << fallbackPort << "...\n";
        closesocket(workerSocket);

        // Try fallback
        workerSocket = socket(AF_INET, SOCK_STREAM, 0);
        workerAddr.sin_port = htons(fallbackPort);

        if (connect(workerSocket, (sockaddr*)&workerAddr, sizeof(workerAddr)) == SOCKET_ERROR)
        {
            cerr << "Server: Fallback worker on port " << fallbackPort << " also not available.\n";
            closesocket(workerSocket);
            return;
        }
    }

    send(workerSocket, inputGene.c_str(), inputGene.size(), 0);

    char buffer[bufferCapacity] = { 0 };
    int bytesReceived = recv(workerSocket, buffer, bufferCapacity - 1, 0);
    if (bytesReceived > 0)
    {
        buffer[bytesReceived] = '\0';
        send(clientSocket, buffer, bytesReceived, 0);
    }

    closesocket(workerSocket);
}

void handleClient(SOCKET clientSocket, int workerPort)
{
    char buffer[bufferCapacity] = { 0 };
    int bytes = recv(clientSocket, buffer, bufferCapacity - 1, 0);
    if (bytes <= 0)
    {
        cerr << "Server: Failed to receive from client.\n";
        closesocket(clientSocket);
        return;
    }

    buffer[bytes] = '\0';
    string inputGene(buffer);
    cout << "Server: Received DNA from client: " << inputGene << endl;

    forwardToWorker(inputGene, clientSocket, workerPort);

    closesocket(clientSocket);
}

void listenPortThread(int listenPort, int workerPort)
{
    SOCKET serverSocket = socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(listenPort);

    bind(serverSocket, (sockaddr*)&addr, sizeof(addr));
    listen(serverSocket, 5);

    cout << "Server: Listening for clients on port " << listenPort << "...\n";

    while (true)
    {
        SOCKET clientSock = accept(serverSocket, nullptr, nullptr);
        if (clientSock != INVALID_SOCKET)
        {
            thread t(handleClient, clientSock, workerPort);
            t.detach();
        }
    }

    closesocket(serverSocket);
}

int main()
{
    WSADATA wsa;
    WSAStartup(MAKEWORD(2, 2), &wsa);

    // Thread for each client port
    thread t1(listenPortThread, 9000, 8000);
    thread t2(listenPortThread, 9001, 8001);

    t1.join();
    t2.join();

    WSACleanup();
    return 0;
}
