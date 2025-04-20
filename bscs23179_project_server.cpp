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

    const char* hGene[][2] =
    {
     {"BRCA1", "ATGAAAAAACTGAGTAAGGAAAGCCTGAGCCAGAGGGTGTCCCGGCTGGGGCATGTGGAGGGTGACTGT"},
     {"TP53", "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGT"},
     {"CFTR", "ATGTTCGTCTTCCTGGATTATGCCTGGCACCATTAAAGAAAATATCATCTTTGGTGTTTCCTATGATGAACACTTGGTTGGC"},
     {"HBB", "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAG"},
     {"INS", "ATGCCCTGTGGATGCGCCTCCTGCACCCCAGGCTTTTGTCAAACAGCACCTTTGTGGTTCTCACTGGTGGGCGCTCAGCCTATCTTG"},
     {"APOE", "ATGAGGCCAGAGGGTCCAGGAGGAAGGTGAGTGAAGAGGGAGTGGAGGGAAGAGGAAGGGAAGGGAGGGAAGAGGAAGGGAGGG"},
     {"BRCA2", "ATGGAGGAGCTCGAGTCGAGGAAGGAGGAGGAGAGGAGGAGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGG"},
     {"EGFR", "ATGAAAGAGGAGGAGGAGGGAGGAGGAAGGAGGGAGGAGGGAGGGAGAGGGAGGGAGGAGGAGGGAGGGAGGGAGGAGGGAGGAAGG"},
     {"FTO", "ATGGAGAGGAGAGGAAGAGGGAGAGGAAGGAGGGAGAGGAAGGGAGGGAGGAGGAGGAGGAAGGAGGAAGGGAGGGAGGGAGGGAGGAGGGAGG"},
     {"MTOR", "ATGAGAGGAGGGAGGGAGGGAGGGAGGGAGGAGGGAGGAGGGAGGGAGGGAGGGAGGGAGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGG"},
     {"VHL", "ATGGAGGAGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGG"},
     {"G6PD", "ATGCCGCAGGGCCTCTCAGGCCAGGGCCTCCTGCGGCTCTGGCCTCTCGGCAGCGTGCAGGCCTCTCTGCCAGGGGCGTCTC"},
     {"DMD", "ATGGATGTCATTTGGGAAAGGAGAGGTGGATGAAGTGAAGGAGGAGGAGGAAGAGGAAGGAGGGAGAGGAGGAGGGGAGGAGGA"},
     {"MYH7", "ATGCCTCGTGCGGCGGTCTCCTGCTGCCTCCTGCCTGCTGCTCCTGCTGCCTGCTCCTGCTGCTGCTGCTGCTGCAGCTG"},
     {"LDLR", "ATGGAGACAGCAGCAGGCGGGGAGGCGGCGCAGCGGGAGGGCGAGGCGCAGCGGCGGAGCGGCGCAGCGGGAGCGGC"},
     {"P53BP1", "ATGGAGGTAGGAGGAAGGAGGAAGGAAGGAAGGAAGGAGGAGGAGGAAGGAGGAAGGAAGGAGGAAGGAGGAGG"},
     {"TTN", "ATGGCCGGAGGTGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAG"},
     {"COL1A1", "ATGGCCTCCCTGGTGGAGCTGGAGGGGCTGGAGGTGGAGGAGGCTGGAGGAGGAGGCTGGAGGAGGAGGAG"},
     {"ACTN3", "ATGGCCGAGCGGCGGAGCGGCGGCGGAGCGGCGGAGCGGCGGCGGAGCGGCGGCGGAGCGGCGGCGGAGC"},
     {"VEGFA", "ATGAACTTTCTGCTGTCTTGGGTGCATTGGAGCCTTGCCTTGCTGCTCTACCTCCACCATGCCAAGTGGTCCCAGGCTGC"}
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
        {"BRCA1", "ACTGAGTAAGGAAAGCCTGAGCCAGAGGGTGTGTGTGTCCAGTTTCTGTTCTTGCAG", "Breast Cancer"},
        {"TP53", "AGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGCCCTTGTCCAGCTCTGGGGAGAGG", "Li-Fraumeni Syndrome"},
        {"CFTR", "ATGTTCGTCTTCCTGGATTATGCCTGGCACCTGCCGTTTTGATGACGCTTCACTG", "Cystic Fibrosis"},
        {"HBB", "CCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTG", "Sickle Cell Anemia"},
        {"HTT", "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG", "Huntington's Disease"},
        {"FBN1", "TGAGGACTTCAGTGAGAGGAGACTTCCAGAGTTGGCTTCCAGAGTCTTCAGTGCAG", "Marfan Syndrome"},
        {"MLH1", "GTGTAGTAGGAGGAGGTTGGAAGTGGATGAGGAGGCTGTTGTAGATGAGTAGGAG", "Hereditary Nonpolyposis Colorectal Cancer"},
        {"PAH", "GGAGTGGGAGTGGTGGTGTGGTGAGGTGGAGTGAGGCTGTGGGAGTGGCTGTGG", "Phenylketonuria"},
        {"DMD", "TGGTCTGGAGGTGTGGTGTGGAGGAGGTGGTGGAGGTGGTGGTGGTGGTGGAGG", "Duchenne Muscular Dystrophy"},
        {"SMN1", "GCTGAGGGTGTGTGTGGTGGAGGAGGTGGGTGGAGGAGGTGTGGAGGTGAGGTG", "Spinal Muscular Atrophy"},
        {"GBA", "GGAGGTGGAGGAGGAGGTGGGTGTGGGTGGAGGTGTGGAGGAGGTGTGGTGGTG", "Gaucher's Disease"},
        {"NPC1", "CTGCTGAGGAGGGTGGAGGAGGTGGAGGGTGGAGGTGGTGGTGGAGGTGGTGGA", "Niemann-Pick Disease Type C"},
        {"MECP2", "AGCTGAGGAGGTGGAGGAGGTGGTGGAGGAGGTGGGTGGAGGAGGTGGAGGTGG", "Rett Syndrome"},
        {"TSC1", "GTGGGTGGAGGTGGAGGTGGAGGAGGTGGAGGTGGTGGAGGTGGGTGGTGGGTG", "Tuberous Sclerosis"},
        {"PKD1", "GGAGGTGGAGGAGGTGGTGGAGGAGGTGGTGGGTGGAGGTGGAGGGTGGAGGTG", "Polycystic Kidney Disease"}
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

string diseasePredictorAndMutation(const string& dbFile, const string& seq)
{
    
    string Bioresult = "Disease Predictor Results:\n";
    sqlite3* db;
    int rc = sqlite3_open(dbFile.c_str(), &db);
    if (rc)
    {
        return Bioresult + "failed to open database:" + string(sqlite3_errmsg(db)) + "\n";
    }

    const char* Diseasequery = R"(
        SELECT geneName, disease FROM MutationGenes
        WHERE INSTR(?, mutationSequence) > 0;
    )";

    sqlite3_stmt* stmt;
    rc = sqlite3_prepare_v2(db, Diseasequery, -1, &stmt, nullptr);
    if (rc != SQLITE_OK)
    {
        return Bioresult + "Predictor couldnt find results: " + string(sqlite3_errmsg(db)) + "\n";
    }

    sqlite3_bind_text(stmt, 1, seq.c_str(), -1, SQLITE_STATIC);

    bool found = false;

    while (sqlite3_step(stmt) == SQLITE_ROW)
    {
        const char* gene = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        const char* disease = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
        Bioresult += "Mutation found in gene " + string(gene) + " Disease: " + string(disease) + "\n";
        found = true;
    }

    if (!found)
    {
        Bioresult += "No mutations found in the client's input sequence.\n";
    }

    sqlite3_finalize(stmt);
    sqlite3_close(db);

    return Bioresult;
}

void comparisonOfGeneSequence(const string& dbFile, const string& inputGene, SOCKET clientSocket)
{
    sqlite3* db;
    sqlite3_stmt* stmt;
    const char* selectSQL = "SELECT name, sequence FROM HealthyGenes";

    if (sqlite3_open(dbFile.c_str(), &db) != SQLITE_OK)
    {
        cerr << "failed to open database.\n";
        return;
    }

    if (sqlite3_prepare_v2(db, selectSQL, -1, &stmt, nullptr) != SQLITE_OK)
    {
        cerr << "failed to select the statement.\n";
        sqlite3_close(db);
        return;
    }

    int bestScore = -1;
    string bestGene;
    AlignmentScore bestAligned;

    while (sqlite3_step(stmt) == SQLITE_ROW)
    {
        string name = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        string Gseq = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));

        AlignmentScore score = smithWaterman(inputGene, Gseq);

        if (score.score > bestScore)
        {
            bestScore = score.score;
            bestGene = name;
            bestAligned = score;
        }
    }

    sqlite3_finalize(stmt);
    sqlite3_close(db);

    string diseaseResults = diseasePredictorAndMutation(dbFile, inputGene);

   

    string output = "Score: " + to_string(bestAligned.score) + "\n"
        + "Aligned Input: " + bestAligned.aligned_input + "\n"
        + "Aligned Gene: " + bestAligned.aligned_healthy + "\n"
        + "Updations: " + to_string(bestAligned.updations)
        + ", Insertions: " + to_string(bestAligned.insertions)
        + ", Deletions: " + to_string(bestAligned.deletions) + "\n\n"
        + diseaseResults;

    send(clientSocket, output.c_str(), output.size(), 0);
    cout << "BioInformatics results sent to client.\n";
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

    comparisonOfGeneSequence("GenesDatabase.db", clientG, clientSocket);

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

    cout << " Server is listening on port : " << port << "\n";

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
    CreateAndInsertDB("GenesDatabase.db");

    WSADATA wsa;
    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0)
    {
        cerr << " WSAStartup failed." << endl;
        return 1;
    }

    thread t1(listenPorThread, 9000);
    thread t2(listenPorThread, 9001);

    t1.join();
    t2.join();

    WSACleanup();
    return 0;
}


