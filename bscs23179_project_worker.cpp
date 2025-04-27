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

    comparisonOfGeneSequence("GeneDatabase.db", clientG, clientSocket);

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

void CreateAndInsertDB(const string& dbPath)
{
    sqlite3* db;
    char* errMsg = nullptr;

    int rc = sqlite3_open(dbPath.c_str(), &db);
    if (rc)
    {
        cerr << "couldnt open database: " << sqlite3_errmsg(db) << endl;
        return;
    }

    //healthy gene database
    const char* createHgeneTableSQL = R"(
        CREATE TABLE IF NOT EXISTS HealthyGenes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT,
            sequence TEXT,
            UNIQUE(name, sequence)
        );
    )";

    rc = sqlite3_exec(db, createHgeneTableSQL, nullptr, nullptr, &errMsg);

    if (rc != SQLITE_OK)
    {
        cerr << "SQL error (HGenes): " << errMsg << endl;
        sqlite3_free(errMsg);
        sqlite3_close(db);
        return;
    }

    const char* hGene[][2] = {
     {"BRCA1", "ATGAAAAAACTGAGTAAGGAAAGCCTGAGCCAGAGGGTGTCCCGGCTGGGGCATGTGGAGGGTGACTGTATGAAAAAACTGAGTAAGGAAAGCCTGAGCCAGAGGGTGTCCCGGCTGGGGCATGTGGAGGGTGACTGTATGAAAAAACTGAGTAAGGAAAGCCTGAGCCAGAGGGTGTCCCGGCTGGGGCATGTGGAGGGTGACTGT"},
     {"TP53", "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGT"},
     {"CFTR", "ATGTTCGTCTTCCTGGATTATGCCTGGCACCATTAAAGAAAATATCATCTTTGGTGTTTCCTATGATGAACACTTGGTTGGCATGTTCGTCTTCCTGGATTATGCCTGGCACCATTAAAGAAAATATCATCTTTGGTGTTTCCTATGATGAACACTTGGTTGGC"},
     {"HBB", "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAG"},
     {"INS", "ATGCCCTGTGGATGCGCCTCCTGCACCCCAGGCTTTTGTCAAACAGCACCTTTGTGGTTCTCACTGGTGGGCGCTCAGCCTATCTTGATGCCCTGTGGATGCGCCTCCTGCACCCCAGGCTTTTGTCAAACAGCACCTTTGTGGTTCTCACTGGTGGGCGCTCAGCCTATCTTG"},
     {"APOE", "ATGAGGCCAGAGGGTCCAGGAGGAAGGTGAGTGAAGAGGGAGTGGAGGGAAGAGGAAGGGAAGGGAGGGAAGAGGAAGGGAGGGATGAGGCCAGAGGGTCCAGGAGGAAGGTGAGTGAAGAGGGAGTGGAGGGAAGAGGAAGGGAAGGGAGGGAAGAGGAAGGGAGGG"},
     {"BRCA2", "ATGGAGGAGCTCGAGTCGAGGAAGGAGGAGGAGAGGAGGAGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGATGGAGGAGCTCGAGTCGAGGAAGGAGGAGGAGAGGAGGAGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGG"},
     {"EGFR", "ATGAAAGAGGAGGAGGAGGGAGGAGGAAGGAGGGAGGAGGGAGGGAGAGGGAGGGAGGAGGAGGGAGGGAGGGAGGAGGGAGGAAGGATGAAAGAGGAGGAGGAGGGAGGAGGAAGGAGGGAGGAGGGAGGGAGAGGGAGGGAGGAGGAGGGAGGGAGGGAGGAGGGAGGAAGG"},
     {"FTO", "ATGGAGAGGAGAGGAAGAGGGAGAGGAAGGAGGGAGAGGAAGGGAGGGAGGAGGAGGAGGAAGGAGGAAGGGAGGGAGGGAGGGAGGAGGGAGGATGGAGAGGAGAGGAAGAGGGAGAGGAAGGAGGGAGAGGAAGGGAGGGAGGAGGAGGAGGAAGGAGGAAGGGAGGGAGGGAGGGAGGAGGGAGG"},
     {"MTOR", "ATGAGAGGAGGGAGGGAGGGAGGGAGGGAGGAGGGAGGAGGGAGGGAGGGAGGGAGGGAGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGATGAGAGGAGGGAGGGAGGGAGGGAGGGAGGAGGGAGGAGGGAGGGAGGGAGGGAGGGAGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGG"},
     {"VHL", "ATGGAGGAGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGATGGAGGAGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGG"},
     {"G6PD", "ATGCCGCAGGGCCTCTCAGGCCAGGGCCTCCTGCGGCTCTGGCCTCTCGGCAGCGTGCAGGCCTCTCTGCCAGGGGCGTCTCATGCCGCAGGGCCTCTCAGGCCAGGGCCTCCTGCGGCTCTGGCCTCTCGGCAGCGTGCAGGCCTCTCTGCCAGGGGCGTCTC"},
     {"DMD", "ATGGATGTCATTTGGGAAAGGAGAGGTGGATGAAGTGAAGGAGGAGGAGGAAGAGGAAGGAGGGAGAGGAGGAGGGGAGGAGGAATGGATGTCATTTGGGAAAGGAGAGGTGGATGAAGTGAAGGAGGAGGAGGAAGAGGAAGGAGGGAGAGGAGGAGGGGAGGAGGA"},
     {"MYH7", "ATGCCTCGTGCGGCGGTCTCCTGCTGCCTCCTGCCTGCTGCTCCTGCTGCCTGCTCCTGCTGCTGCTGCTGCTGCAGCTGATGCCTCGTGCGGCGGTCTCCTGCTGCCTCCTGCCTGCTGCTCCTGCTGCCTGCTCCTGCTGCTGCTGCTGCTGCAGCTG"},
     {"LDLR", "ATGGAGACAGCAGCAGGCGGGGAGGCGGCGCAGCGGGAGGGCGAGGCGCAGCGGCGGAGCGGCGCAGCGGGAGCGGCATGGAGACAGCAGCAGGCGGGGAGGCGGCGCAGCGGGAGGGCGAGGCGCAGCGGCGGAGCGGCGCAGCGGGAGCGGC"},
     {"P53BP1", "ATGGAGGTAGGAGGAAGGAGGAAGGAAGGAAGGAAGGAGGAGGAGGAAGGAGGAAGGAAGGAGGAAGGAGGAGGATGGAGGTAGGAGGAAGGAGGAAGGAAGGAAGGAAGGAGGAGGAGGAAGGAGGAAGGAAGGAGGAAGGAGGAGG"},
     {"TTN", "ATGGCCGGAGGTGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGATGGCCGGAGGTGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGG"},
     {"COL1A1", "ATGGCCTCCCTGGTGGAGCTGGAGGGGCTGGAGGTGGAGGAGGCTGGAGGAGGAGGCTGGAGGAGGAGGAGATGGCCTCCCTGGTGGAGCTGGAGGGGCTGGAGGTGGAGGAGGCTGGAGGAGGAGGCTGGAGGAGGAGGAG"},
     {"ACTN3", "ATGGCCGAGCGGCGGAGCGGCGGCGGAGCGGCGGAGCGGCGGCGGAGCGGCGGCGGAGCGGCGGCGGAGCATGGCCGAGCGGCGGAGCGGCGGCGGAGCGGCGGAGCGGCGGCGGAGCGGCGGCGGAGCGGCGGCGGAGC"},
     {"VEGFA", "ATGAACTTTCTGCTGTCTTGGGTGCATTGGAGCCTTGCCTTGCTGCTCTACCTCCACCATGCCAAGTGGTCCCAGGCTGCATGAACTTTCTGCTGTCTTGGGTGCATTGGAGCCTTGCCTTGCTGCTCTACCTCCACCATGCCAAGTGGTCCCAGGCTGC"}
    };




    const char* insertHgeneSQL = R"(
    INSERT OR IGNORE INTO HealthyGenes (name, sequence)
    VALUES (?, ?);
    )";

    sqlite3_stmt* insertHgeneStmt;

    rc = sqlite3_prepare_v2(db, insertHgeneSQL, -1, &insertHgeneStmt, nullptr);
    if (rc != SQLITE_OK)
    {
        cerr << "preparation failed at Hgene: " << sqlite3_errmsg(db) << "\n";
        sqlite3_close(db);
        return;
    }

    for (int i = 0; i < sizeof(hGene) / sizeof(hGene[0]); ++i)
    {
        sqlite3_reset(insertHgeneStmt);
        sqlite3_clear_bindings(insertHgeneStmt);
        sqlite3_bind_text(insertHgeneStmt, 1, hGene[i][0], -1, SQLITE_STATIC);
        sqlite3_bind_text(insertHgeneStmt, 2, hGene[i][1], -1, SQLITE_STATIC);

        rc = sqlite3_step(insertHgeneStmt);
        if (rc != SQLITE_DONE)
            cerr << "Insertion failed at Hgene: " << sqlite3_errmsg(db) << endl;
    }

    sqlite3_finalize(insertHgeneStmt);



    // mutation gene database
    const char* createMgeneTableSQL = R"(
        CREATE TABLE IF NOT EXISTS MutationGenes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            geneName TEXT,
            mutationSequence TEXT,
            disease TEXT,
            UNIQUE(geneName, mutationSequence)
        );
    )";

    rc = sqlite3_exec(db, createMgeneTableSQL, nullptr, nullptr, &errMsg);
    if (rc != SQLITE_OK)
    {
        cerr << "SQL error (MGenes): " << errMsg << endl;
        sqlite3_free(errMsg);
        sqlite3_close(db);
        return;
    }


    const char* insertMgeneSQL = R"(
        INSERT OR IGNORE INTO MutationGenes (geneName, mutationSequence, disease)
        VALUES (?, ?, ?);
    )";

    const char* mGene[][3] =
    {
        {"BRCA1", "ACTGAGTAAGGAAAGCCTGAGCCAGAGGGTGTGTGTGTCCAGTTTCTGTTCTTGCAGACTGAGTAAGGAAAGCCTGAGCCAGAGGGTGTGTGTGTCCAGTTTCTGTTCTTGCAGACTGAGTAAGGAAAGCCTGAGCCAGAGGGTGTGTGTGTCCAGTTTCTGTTCTTGCAG", "Breast Cancer"},
        {"TP53", "AGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGCCCTTGTCCAGCTCTGGGGAGAGGAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGCCCTTGTCCAGCTCTGGGGAGAGGAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGCCCTTGTCCAGCTCTGGGGAGAGG", "Li-Fraumeni Syndrome"},
        {"CFTR", "ATGTTCGTCTTCCTGGATTATGCCTGGCACCTGCCGTTTTGATGACGCTTCACTGATGTTCGTCTTCCTGGATTATGCCTGGCACCTGCCGTTTTGATGACGCTTCACTGATGTTCGTCTTCCTGGATTATGCCTGGCACCTGCCGTTTTGATGACGCTTCACTG", "Cystic Fibrosis"},
        {"HBB", "CCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGCCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGCCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTG", "Sickle Cell Anemia"},
        {"HTT", "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG", "Huntington's Disease"},
        {"FBN1", "TGAGGACTTCAGTGAGAGGAGACTTCCAGAGTTGGCTTCCAGAGTCTTCAGTGCAGTGAGGACTTCAGTGAGAGGAGACTTCCAGAGTTGGCTTCCAGAGTCTTCAGTGCAGTGAGGACTTCAGTGAGAGGAGACTTCCAGAGTTGGCTTCCAGAGTCTTCAGTGCAG", "Marfan Syndrome"},
        {"MLH1", "GTGTAGTAGGAGGAGGTTGGAAGTGGATGAGGAGGCTGTTGTAGATGAGTAGGAGGTGTAGTAGGAGGAGGTTGGAAGTGGATGAGGAGGCTGTTGTAGATGAGTAGGAGGTGTAGTAGGAGGAGGTTGGAAGTGGATGAGGAGGCTGTTGTAGATGAGTAGGAG", "Hereditary Nonpolyposis Colorectal Cancer"},
        {"PAH", "GGAGTGGGAGTGGTGGTGTGGTGAGGTGGAGTGAGGCTGTGGGAGTGGCTGTGGGGAGTGGGAGTGGTGGTGTGGTGAGGTGGAGTGAGGCTGTGGGAGTGGCTGTGGGGAGTGGGAGTGGTGGTGTGGTGAGGTGGAGTGAGGCTGTGGGAGTGGCTGTGG", "Phenylketonuria"},
        {"DMD", "TGGTCTGGAGGTGTGGTGTGGAGGAGGTGGTGGAGGTGGTGGTGGTGGTGGAGGTGGTCTGGAGGTGTGGTGTGGAGGAGGTGGTGGAGGTGGTGGTGGTGGTGGAGGTGGTCTGGAGGTGTGGTGTGGAGGAGGTGGTGGAGGTGGTGGTGGTGGTGGAGG", "Duchenne Muscular Dystrophy"},
        {"SMN1", "GCTGAGGGTGTGTGTGGTGGAGGAGGTGGGTGGAGGAGGTGTGGAGGTGAGGTGGCTGAGGGTGTGTGTGGTGGAGGAGGTGGGTGGAGGAGGTGTGGAGGTGAGGTGGCTGAGGGTGTGTGTGGTGGAGGAGGTGGGTGGAGGAGGTGTGGAGGTGAGGTG", "Spinal Muscular Atrophy"},
        {"GBA", "GGAGGTGGAGGAGGAGGTGGGTGTGGGTGGAGGTGTGGAGGAGGTGTGGTGGTGGGAGGTGGAGGAGGAGGTGGGTGTGGGTGGAGGTGTGGAGGAGGTGTGGTGGTGGGAGGTGGAGGAGGAGGTGGGTGTGGGTGGAGGTGTGGAGGAGGTGTGGTGGTG", "Gaucher's Disease"},
        {"NPC1", "CTGCTGAGGAGGGTGGAGGAGGTGGAGGGTGGAGGTGGTGGTGGAGGTGGTGGACTGCTGAGGAGGGTGGAGGAGGTGGAGGGTGGAGGTGGTGGTGGAGGTGGTGGACTGCTGAGGAGGGTGGAGGAGGTGGAGGGTGGAGGTGGTGGTGGAGGTGGTGGA", "Niemann-Pick Disease Type C"},
        {"MECP2", "AGCTGAGGAGGTGGAGGAGGTGGTGGAGGAGGTGGGTGGAGGAGGTGGAGGTGGAGCTGAGGAGGTGGAGGAGGTGGTGGAGGAGGTGGGTGGAGGAGGTGGAGGTGGAGCTGAGGAGGTGGAGGAGGTGGTGGAGGAGGTGGGTGGAGGAGGTGGAGGTGG", "Rett Syndrome"},
        {"TSC1", "GTGGGTGGAGGTGGAGGTGGAGGAGGTGGAGGTGGTGGAGGTGGGTGGTGGGTGGTGGGTGGAGGTGGAGGTGGAGGAGGTGGAGGTGGTGGAGGTGGGTGGTGGGTGGTGGGTGGAGGTGGAGGTGGAGGAGGTGGAGGTGGTGGAGGTGGGTGGTGGGTG", "Tuberous Sclerosis"},
        {"PKD1", "GGAGGTGGAGGAGGTGGTGGAGGAGGTGGTGGGTGGAGGTGGAGGGTGGAGGTGGGAGGTGGAGGAGGTGGTGGAGGAGGTGGTGGGTGGAGGTGGAGGGTGGAGGTGGGAGGTGGAGGAGGTGGTGGAGGAGGTGGTGGGTGGAGGTGGAGGGTGGAGGTG", "Polycystic Kidney Disease"},
        {"ALB", "ATGGCTGAGAACAGTCACAGTGTGAGGCTGTGGTTGCCTTGGGTGTGTGGGCCTGGTGCTGCTGAGGCTGGGAGGAGGGGCAGGAGG", "Analbuminemia"},
        {"TTR", "ATGGAAGTGTTTGGGTTGGGGTTGCTGAGGTGGGAGTGGGGCTTCTGGTGTGGGGTTGGTGCTGTGGTGGTTGGTGGTTTGGGT", "Familial Amyloid Polyneuropathy"},
        {"HFE", "ATGGTGCTGTGGTGGGGCTGCTGAGGAGGAGGTGAGGGAGGGGGAGGTGGGGGTGAGGGAGGGGGAGGTGGGGGTG", "Hereditary Hemochromatosis"},
        {"COL6A1", "ATGCCAGGGTGGGGAGGGTGGGTGAGGGAGGGAGGAGGAGGGAGGAGGGGAGGGTGGAGGAGGGAGGGTGGAGGG", "Bethlem Myopathy"},
        {"LMNA", "ATGGAGACAGAGGAGGAGGGAGAGGAGGGAGGGAGGAGGGAGGAGGGAGGAGGGAGGGAGGAGGGAGGGAGGA", "Emery–Dreifuss Muscular Dystrophy"},
        {"FGFR3", "ATGCGCAGCTTCCCGGAGAGGAGAGAGAGAGAGCTGGTGGAGGAGGAGGGAGGGAGGGAGGGAGGGAGGAGGA", "Achondroplasia"},
        {"PCSK9", "ATGGAGGGAGGCGGAGGGGAGGGGAGGGAGGGGAGGGAGGGGAGGGGAGGGAGGGAGGGAGGGGAGGGAGGGA", "Hypercholesterolemia"}
    };



    sqlite3_stmt* insertMgeneStmt;
    rc = sqlite3_prepare_v2(db, insertMgeneSQL, -1, &insertMgeneStmt, nullptr);
    if (rc != SQLITE_OK)
    {
        cerr << "preparation failed Mgene: " << sqlite3_errmsg(db) << endl;
        sqlite3_close(db);
        return;
    }

    for (int i = 0; i < sizeof(mGene) / sizeof(mGene[0]); ++i)
    {
        sqlite3_reset(insertMgeneStmt);
        sqlite3_clear_bindings(insertMgeneStmt);
        sqlite3_bind_text(insertMgeneStmt, 1, mGene[i][0], -1, SQLITE_STATIC);
        sqlite3_bind_text(insertMgeneStmt, 2, mGene[i][1], -1, SQLITE_STATIC);
        sqlite3_bind_text(insertMgeneStmt, 3, mGene[i][2], -1, SQLITE_STATIC);

        rc = sqlite3_step(insertMgeneStmt);
        if (rc != SQLITE_DONE)
            cerr << "Insertion failed at Mgene : " << sqlite3_errmsg(db) << endl;
    }

    sqlite3_finalize(insertMgeneStmt);
    sqlite3_close(db);
}




void handleServerRequest(SOCKET serverSocket)
{
    char buffer[bufferCapacity] = { 0 };
    int recvBytes = recv(serverSocket, buffer, bufferCapacity - 1, 0);
    if (recvBytes <= 0)
    {
        cerr << "Worker: Failed to receive sequence.\n";
        closesocket(serverSocket);
        return;
    }

    buffer[recvBytes] = '\0';
    string inputGene(buffer);

    cout << "Worker: Received gene from server: " << inputGene << endl;

    sqlite3_initialize();
    sqlite3_shutdown(); // Just safety; actual DB call already opens/uses sqlite3

    SOCKET dummyClient = serverSocket; // Reuse socket for response
    comparisonOfGeneSequence("GeneDatabase.db", inputGene, dummyClient);

    closesocket(serverSocket);
}

int main()
{
    CreateAndInsertDB("GeneDatabase.db");
    WSADATA wsa;
    WSAStartup(MAKEWORD(2, 2), &wsa);

    SOCKET workerSocket = socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in workerAddr{};
    workerAddr.sin_family = AF_INET;
    workerAddr.sin_addr.s_addr = INADDR_ANY;
    workerAddr.sin_port = htons(8000); // Worker port

    bind(workerSocket, (sockaddr*)&workerAddr, sizeof(workerAddr));
    listen(workerSocket, 5);

    cout << "Worker: Listening for DNA tasks from server...\n";

    while (true)
    {
        SOCKET serverSock = accept(workerSocket, nullptr, nullptr);
        if (serverSock != INVALID_SOCKET)
        {
            handleServerRequest(serverSock);
        }
    }

    closesocket(workerSocket);
    WSACleanup();
    return 0;
}




