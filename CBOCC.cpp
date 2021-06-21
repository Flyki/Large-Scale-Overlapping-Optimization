#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <sstream>
#include "Header.h"
#include <float.h>
#include <cstdlib>
#include "CBOG_CBD.h"

using namespace std;

Benchmarks* generateFuncObj(int funcID) {
	Benchmarks* fp;
	if (funcID == 1) {
		fp = new F1();
	}
	else if (funcID == 2) {
		fp = new F2();
	}
	else if (funcID == 3) {
		fp = new F3();
	}
	else if (funcID == 4) {
		fp = new F4();
	}
	else if (funcID == 5) {
		fp = new F5();
	}
	else if (funcID == 6) {
		fp = new F6();
	}
	else if (funcID == 7) {
		fp = new F7();
	}
	else if (funcID == 8) {
		fp = new F8();
	}
	else if (funcID == 9) {
		fp = new F9();
	}
	else if (funcID == 10) {
		fp = new F10();
	}
	else if (funcID == 11) {
		fp = new F11();
	}
	else if (funcID == 12) {
		fp = new F12();
	}
	else {
		cerr << "Fail to locate Specified Function Index" << endl;
		exit(-1);
	}
	return fp;
}

int main(int argc, char* argv[]) {
	int func = atoi(argv[1]);
	string method(argv[2]);
	int seed = atoi(argv[3]);
	long int maxfes = atol(argv[4]);
	int taskId = seed;
	int DIM = 905;
	Benchmarks* fp;
	fp = generateFuncObj(func);
	ifstream groupfile;
	vector<vector<int>> groups;
	stringstream ss;
	string funcStr;
	ss << func;
	ss >> funcStr;
	ss.clear();
	string fileName = funcStr + "po.txt";
	groupfile.open(fileName.c_str());
	for (int i = 0; i < 20; i++) {
		int elenum;
		vector<int> v;
		groupfile >> elenum;
		for (int j = 0; j < elenum; j++) {
			int ele;
			groupfile >> ele;
			v.push_back(ele);
		}
		groups.push_back(v);
	}
	groupfile.close();
	fileName = funcStr + "oo.txt";
	ifstream overfile(fileName.c_str());
	vector<vector<int>> overiables;
	vector<vector<int>> overiablesRedandunt;
	map<int, vector<pair<int, int>>> sharedvar_group_pos;
	set<int> allsharedv;
	for (int i = 0; i < 20; i++) {
		int eleNum;
		vector<int> v;
		vector<int> v2;
		overfile >> eleNum;
		for (int j = 0; j < eleNum; j++) {
			int ele;
			overfile >> ele;
			if (allsharedv.find(ele) == allsharedv.end()) {
				v.push_back(ele);
				allsharedv.insert(ele);
			}
			v2.push_back(ele);
		}
		if (!v.empty()) {
			overiables.push_back(v);
		}
		overiablesRedandunt.push_back(v2);
	}
	overfile.close();
	for (int groupn = 0; groupn < (int)groups.size(); groupn++) {
		for (int i = 0; i < (int)groups[groupn].size(); i++) {
			for (int j = 0; j < (int)overiablesRedandunt[groupn].size(); j++) {
				if (groups[groupn][i] == overiablesRedandunt[groupn][j]) {
					sharedvar_group_pos[overiablesRedandunt[groupn][j]].push_back(make_pair(groupn, i));
				}
			}
		}
	}

	if (method == "CBCCO") {
		CBOG_CBD oneSolver(fp, DIM, seed, 100, groups, overiables, overiablesRedandunt, sharedvar_group_pos, maxfes);
		oneSolver.testStage();
		oneSolver.optimizationStage();
	}
	else {
		cout << "No such way!!" << endl;
	}
	delete fp;
	return 0;
}