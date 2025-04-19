#include <iostream>
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

};


AlignmentScore smithWaterman( const string& input, const string& healthy, int match = 2,int mismatch = -1, int gap = -2 )
{
    int o = input.size();
    int h = healthy.size();

    vector<vector<int>> dp(o + 1, vector<int>(h + 1, 0));
    vector<vector<int>> traceback(o + 1, vector<int>(h + 1, 0));

    int max1 = 0;
    int max2 = 0;
    int maxScore = 0;
    

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
                traceback[i][j] = 1;
            }
            else if (dp[i][j] == UpScore)
            {
                traceback[i][j] = 2;
            }
            else if (dp[i][j] == LeftScore)
            {
                traceback[i][j] = 3;
            }

            if (dp[i][j] > maxScore) 
            {
                maxScore = dp[i][j];
                max1 = i;
                max2 = j;
            }
        }
    }
    
  

    return { maxScore };
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

    
    AlignmentScore score = smithWaterman(clientG, healthyG);

 
    string output = "Alignment Score: " + std::to_string(score.score) + "\n";
    send(clientSocket, output.c_str(), output.size(), 0);
    cout << "Alignment score send to client: " << score.score << endl;

    closesocket(clientSocket);
}

void listenPorThread(int port)
{
    SOCKET ss = socket(AF_INET, SOCK_STREAM, 0);
    if (ss == INVALID_SOCKET)
    {
        cerr << " Socket couldnot built."<<endl;
        return;
    }

    sockaddr_in server{};
    server.sin_family = AF_INET;
    server.sin_addr.s_addr = INADDR_ANY;
    server.sin_port = htons(port);

    if (bind(ss, (sockaddr*)&server, sizeof(server)) == SOCKET_ERROR)
    {
        cerr << " Bind failed."<< endl;
        closesocket(ss);
        return;
    }

    if (listen(ss, 5) == SOCKET_ERROR)
    {
        cerr << " Listening failed.\n";
        closesocket(ss);
        return;
    }

    cout << " Server is listening on port " << port << endl;

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
