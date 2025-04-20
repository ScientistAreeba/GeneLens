#include <iostream>
#include "sqlite3.h"
#define NOMINMAX
#include <winsock2.h>
#include <ws2tcpip.h>
#include <thread>
#include <string>
#include <vector>
using namespace std;

#pragma comment(lib, "ws2_32.lib")
#define bufferCapacity 1024


string healthyG = "ATGCGTAGTC";

struct AlignmentScore
{
    int score;
    string aligned_input;
    string aligned_healthy;
    int updations;
    int insertions;
    int deletions;
};


AlignmentScore smithWaterman(const string& input, const string& healthy, int match = 2, int mismatch = -1, int gap = -2)
{
    int o = input.size();
    int h = healthy.size();
    vector<vector<int>> dp(o + 1, vector<int>(h + 1, 0));
    vector<vector<int>> tracing(o + 1, vector<int>(h + 1, 0));

    int maxScore = 0;
    int max1 = 0;
    int max2 = 0;


    for (int i = 1; i <= o; ++i)
    {
        for (int j = 1; j <= h; ++j)
        {

            int diagonalScore;
            if (input[i - 1] == healthy[j - 1])
            {
                diagonalScore = dp[i - 1][j - 1] + match;
            }
            else {
                diagonalScore = dp[i - 1][j - 1] + mismatch;
            }


            int UpScore = dp[i - 1][j] + gap;
            int LeftScore = dp[i][j - 1] + gap;

            dp[i][j] = max({ 0, diagonalScore, UpScore, LeftScore });

            if (dp[i][j] == diagonalScore)
            {
                tracing[i][j] = 1;
            }
            else if (dp[i][j] == UpScore)
            {
                tracing[i][j] = 2;
            }
            else if (dp[i][j] == LeftScore)
            {
                tracing[i][j] = 3;
            }

            if (dp[i][j] > maxScore)
            {
                maxScore = dp[i][j];
                max1 = i;
                max2 = j;
            }
        }
    }


    string alignInput;
    string alignHealthy;
    int i = max1;
    int j = max2;

    int update = 0;
    int insert = 0;
    int remove = 0;


    while (i > 0 && j > 0 && dp[i][j] != 0)
    {
        if (tracing[i][j] == 1)
        {
            alignInput = input[i - 1] + alignInput;
            alignHealthy = healthy[j - 1] + alignHealthy;
            if (input[i - 1] != healthy[j - 1]) ++update;
            --i; --j;
        }
        else if (tracing[i][j] == 2)
        {
            alignInput = input[i - 1] + alignInput;
            alignHealthy = "-" + alignHealthy;
            ++insert;
            --i;
        }
        else if (tracing[i][j] == 3)
        {
            alignInput = "-" + alignInput;
            alignHealthy = healthy[j - 1] + alignHealthy;
            ++remove;
            --j;
        }
        else
        {
            break;
        }
    }

    return { maxScore, alignInput, alignHealthy, update, insert, remove };

}



void clientComputation(SOCKET clientSocket)
{
    char bufferSpace[bufferCapacity] = { 0 };
    int inputSize = recv(clientSocket, bufferSpace, bufferCapacity - 1, 0);

    if (inputSize <= 0)
    {
        cerr << "Failed to receive gene sequence.\n";
        closesocket(clientSocket);
        return;
    }

    bufferSpace[inputSize] = '\0';
    string clientG(bufferSpace);

    cout << "\n Received DNA from client: " << clientG << "\n";

    //compareWithHealthyGenesFromDB("GenesDatabase.db", clientG, clientSocket);

    closesocket(clientSocket);
}

void listenPorThread(int port)
{
    SOCKET ss = socket(AF_INET, SOCK_STREAM, 0);
    if (ss == INVALID_SOCKET)
    {
        cerr << " Socket couldnot built." << endl;
        return;
    }

    sockaddr_in server{};
    server.sin_family = AF_INET;
    server.sin_addr.s_addr = INADDR_ANY;
    server.sin_port = htons(port);

    if (bind(ss, (sockaddr*)&server, sizeof(server)) == SOCKET_ERROR)
    {
        cerr << " Bind failed." << endl;
        closesocket(ss);
        return;
    }

    if (listen(ss, 5) == SOCKET_ERROR)
    {
        cerr << " Listening failed.\n";
        closesocket(ss);
        return;
    }

    cout << " Server is listening on port \n" << port << endl;

    while (true)
    {
        sockaddr_in client{};
        int client_size = sizeof(client);
        SOCKET client_socket = accept(ss, (sockaddr*)&client, &client_size);

        if (client_socket != INVALID_SOCKET)
        {
            thread client_thread(clientComputation, client_socket);
            client_thread.detach();
        }
    }

    closesocket(ss);
}


int main()
{
    WSADATA wsa;


    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0)
    {
        cerr << " WSAStartup failed."<< endl;
        return 1;
    }

    thread t1(listenPorThread, 9000);
    thread t2(listenPorThread, 9001);

    t1.join();
    t2.join();

    WSACleanup();
    return 0;
}
