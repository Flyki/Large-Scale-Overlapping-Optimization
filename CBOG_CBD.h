#pragma once
#include "CMAESO.h"
#include <cfloat>
#include <random>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <set>
#include <ctime>
#include "Header.h"

class CBOG_CBD
{
public:
	CBOG_CBD(Benchmarks*, int DIM, int seed, int testRound, vector<vector<int>> groups, vector<vector<int>> overiables, vector<vector<int>> overiablesRedundant, map<int, vector<pair<int, int>>> share_group_pos, long int MAXFES);
	~CBOG_CBD();
	void initializeOptimizers();
	void testStage();
	void storeInformation(int);
	void writeRestartFile(vector<int>, int);
	void writeRestartFile2(vector<int>, int);
	void optimizationStage();
	vector<int> getAwardList();

	Benchmarks* fp;
	vector<CMAESO*> optimizers;
	double* contribution;
	double* gbest;
	double gbestf;
	long int usedFEs;
	long int MAXFES;

	double** allC;
	double* allmean;
	double* allps;
	double* allpc;
	double* allsigma;

	int DIM;
	int seed;
	int testRound;

	vector<vector<int>> groups;
	vector<vector<int>> overiables;
	vector<vector<int>> overiablesRedundant;
	map<int, vector<pair<int, int>>> share_group_pos;
	ofstream result;

	clock_t starttime;
	clock_t endtime;
};
