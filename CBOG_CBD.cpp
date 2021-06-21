#include "CBOG_CBD.h"

CBOG_CBD::CBOG_CBD(Benchmarks* fp, int DIM, int seed, int testRound, vector<vector<int>> groups, vector<vector<int>> overiables, vector<vector<int>> overiablesRedundant, map<int, vector<pair<int, int>>> share_group_pos, long int MAXFES) {
	starttime = clock();
	this->groups = groups;
	this->overiables = overiables;
	this->overiablesRedundant = overiablesRedundant;
	this->share_group_pos = share_group_pos;
	this->fp = fp;
	this->DIM = DIM;
	this->seed = seed;
	this->usedFEs = 0;
	this->MAXFES = MAXFES;
	this->testRound = testRound;
	this->gbest = new double[this->DIM];
	for (int i = 0; i < this->DIM; i++) {
		this->gbest[i] = ((double)this->fp->getMaxX() + this->fp->getMinX()) / 2.0;
	}
	this->gbestf = fp->compute(gbest);
	this->allmean = new double[DIM];
	for (int i = 0; i < this->DIM; i++) {
		allmean[i] = ((double)this->fp->getMaxX() + this->fp->getMinX()) / 2.0;
	}
	this->allC = new double* [DIM];
	for (int i = 0; i < DIM; i++) {
		this->allC[i] = new double[DIM];
		memset(allC[i], 0, sizeof(double) * DIM);
		allC[i][i] = 1;
	}
	this->allps = new double[DIM];
	this->allpc = new double[DIM];
	memset(allps, 0, sizeof(double) * DIM);
	memset(allpc, 0, sizeof(double) * DIM);
	int groupn = groups.size();
	this->allsigma = new double[groupn];
	for (int i = 0; i < groupn; i++) {
		allsigma[i] = ((double)this->fp->getMaxX() - this->fp->getMinX()) / 3.0;
	}
	this->contribution = new double[groupn];
	memset(contribution, 0, sizeof(double) * groupn);
	stringstream ss;
	ss << fp->getID();
	string resultfilename, funcStr, seedStr;
	ss >> funcStr;
	ss.clear();
	ss << seed;
	ss >> seedStr;
	ss.clear();
	ss << testRound;
	string oneRstr;
	ss >> oneRstr;
	resultfilename = funcStr + "." + seedStr + "." + oneRstr + ".CBOG-CBD.result.txt";
	result.open(resultfilename.c_str(), ios::app);

	initializeOptimizers();
}

void CBOG_CBD::initializeOptimizers() {
	for (int i = 0; i < (int)groups.size(); i++) {
		vector<int> dimtemp = groups[i];
		vector<int>::iterator itr = dimtemp.begin();
		while (itr != dimtemp.end())
		{
			if (share_group_pos.find(*itr) == share_group_pos.end()) {
				itr++;
			}
			else {
				dimtemp.erase(itr);
				itr = dimtemp.begin();
			}
		}
		CMAESO* oneoptimizer = new CMAESO(fp, i, dimtemp, dimtemp.size() * 2, seed);
		this->optimizers.push_back(oneoptimizer);
	}
}

void CBOG_CBD::testStage() {
	memset(contribution, 0, sizeof(double) * groups.size());
	for (int rr = 0; rr < testRound; rr++) {
		for (int i = 0; i < (int)optimizers.size(); i++) {
			double oldf = gbestf;
			optimizers[i]->generateSubsolution();
			optimizers[i]->calculateFitness(gbest, DIM, usedFEs);
			optimizers[i]->updateTheDistribution();
			optimizers[i]->updateContextVector(gbest, gbestf, DIM);
			contribution[i] += (oldf - gbestf);
		}
		result << usedFEs << "," << gbestf << endl;
	}
	for (int i = 0; i < (int)optimizers.size(); i++) {
		storeInformation(i);
	}
}

void CBOG_CBD::optimizationStage() {
	int* overlapAssign = new int[overiables.size()];
	vector<vector<int>> tempgroups = groups;
	for (int i = 0; i < (int)overiables.size(); i++) {
		int theOneGonnaErase = i;
		if (contribution[i] > contribution[i + 1])
			theOneGonnaErase = i + 1;
		overlapAssign[i] = i + i + 1 - theOneGonnaErase;
		for (int j = 0; j < (int)overiables[i].size(); j++) {
			vector<int>::iterator itr = tempgroups[theOneGonnaErase].begin();
			while (itr != tempgroups[theOneGonnaErase].end())
			{
				if (*itr == overiables[i][j]) {
					tempgroups[theOneGonnaErase].erase(itr);
					break;
				}
				itr++;
			}
		}
	}
	for (int i = 0; i < (int)optimizers.size(); i++) {
		writeRestartFile(tempgroups[i], i);
		optimizers[i]->resume(tempgroups[i], 0);
	}
	int rr = 0;
	double* oldgbest = new double[DIM];
	memcpy(oldgbest, gbest, sizeof(double) * DIM);
	while (usedFEs < MAXFES)
	{
		rr++;
		//round-robin optimization
		for (int i = 0; i < (int)optimizers.size(); i++) {
			double oldf = gbestf;
			optimizers[i]->generateSubsolution();
			optimizers[i]->calculateFitness(gbest, DIM, usedFEs);
			optimizers[i]->updateTheDistribution();
			optimizers[i]->updateContextVector(gbest, gbestf, DIM);
			contribution[i] = 0.5 * contribution[i] + 0.5 * (oldf - gbestf);
		}

		//award some subcomponents
		vector<int> awardlist = getAwardList();
		for (auto const e : awardlist) {
			double oldf = gbestf;
			optimizers[e]->generateSubsolution();
			optimizers[e]->calculateFitness(gbest, DIM, usedFEs);
			optimizers[e]->updateTheDistribution();
			optimizers[e]->updateContextVector(gbest, gbestf, DIM);
			contribution[e] = 0.5 * contribution[e] + 0.5 * (oldf - gbestf);
		}

		result << usedFEs << "," << gbestf << endl;
	}
	delete[] overlapAssign;
	delete[] oldgbest;
}

vector<int> CBOG_CBD::getAwardList() {
	vector<int> contrindexRank;
	for (int i = 0; i < (int)optimizers.size(); i++) {
		contrindexRank.push_back(i);
	}
	for (int i = 0; i < (int)optimizers.size() - 1; i++) {
		for (int j = 0; j < (int)optimizers.size() - 1 - i; j++) {
			if (contribution[contrindexRank[j]] < contribution[contrindexRank[j + 1]]) {
				int temp = contrindexRank[j];
				contrindexRank[j] = contrindexRank[j + 1];
				contrindexRank[j + 1] = temp;
			}
		}
	}
	vector<int> awardlist;
	double biggestContri = contribution[contrindexRank[0]];
	for (int i = 0; i < (int)optimizers.size(); i++) {
		if (biggestContri / contribution[contrindexRank[i]] < 2.0) {
			awardlist.push_back(contrindexRank[i]);
		}
		else {
			break;
		}
	}
	if (awardlist.size() == optimizers.size()) {
		awardlist.clear();
	}
	sort(awardlist.begin(), awardlist.end());
	return awardlist;
}

void CBOG_CBD::storeInformation(int index) {
	//store mean
	//const double* xmean = cmaes_GetPtr(&(optimizers[index]->evo), "xmean");
	for (int i = 0; i < optimizers[index]->novar; i++) {
		allmean[optimizers[index]->dims[i]] = optimizers[index]->evo.rgxmean[i];
	}
	//store ps
	for (int i = 0; i < optimizers[index]->novar; i++) {
		allps[optimizers[index]->dims[i]] = optimizers[index]->evo.rgps[i];
	}
	//store pc
	for (int i = 0; i < optimizers[index]->novar; i++) {
		allpc[optimizers[index]->dims[i]] = optimizers[index]->evo.rgpc[i];
	}
	//store C
	for (int i = 0; i < optimizers[index]->novar; i++) {
		for (int j = 0; j <= i; j++) {
			allC[optimizers[index]->dims[i]][optimizers[index]->dims[j]] = optimizers[index]->evo.C[i][j];
			allC[optimizers[index]->dims[j]][optimizers[index]->dims[i]] = optimizers[index]->evo.C[i][j];
		}
	}
	//store sigma
	allsigma[index] = optimizers[index]->evo.sigma;
}

void CBOG_CBD::writeRestartFile(vector<int> temps, int index) {
	stringstream ss;
	ss << fp->getID();
	string funcStr;
	ss >> funcStr;
	ss.clear();
	ss << index;
	string idstr;
	ss >> idstr;
	string fileName = funcStr + "_" + idstr + "_resume.par";
	ofstream statefile(fileName.c_str());
	statefile << '\n';
	statefile << "# resume " << temps.size() << '\n';
	statefile << "xmean" << '\n';
	for (int i = 0; i < (int)temps.size(); i++) {
		statefile << allmean[temps[i]] << (i == (temps.size() - 1) ? "\n" : "\t");
	}
	statefile << "path for sigma\n";
	for (int i = 0; i < (int)temps.size(); i++) {
		statefile << allps[temps[i]] << (i == (temps.size() - 1) ? "\n" : "\t");
	}
	statefile << "path for C\n";
	for (int i = 0; i < (int)temps.size(); i++) {
		statefile << allpc[temps[i]] << (i == (temps.size() - 1) ? "\n" : "\t");
	}
	statefile << "sigma " << ((double)fp->getMaxX() - fp->getMinX()) / 3.0 << '\n';
	statefile << "covariance matrix\n";
	for (int i = 0; i < (int)temps.size(); i++) {
		for (int j = 0; j <= i; j++) {
			statefile << allC[temps[i]][temps[j]] << (j == i ? "\n" : "\t");
		}
	}
	statefile.close();
}

void CBOG_CBD::writeRestartFile2(vector<int> temps, int index) {
	stringstream ss;
	ss << fp->getID();
	string funcStr;
	ss >> funcStr;
	ss.clear();
	ss << index;
	string idstr;
	ss >> idstr;
	string fileName = funcStr + "_" + idstr + "_resume.par";
	ofstream statefile(fileName.c_str());
	statefile << '\n';
	statefile << "# resume " << temps.size() << '\n';
	statefile << "xmean" << '\n';
	for (int i = 0; i < (int)temps.size(); i++) {
		statefile << allmean[temps[i]] << (i == (temps.size() - 1) ? "\n" : "\t");
	}
	statefile << "path for sigma\n";
	for (int i = 0; i < (int)temps.size(); i++) {
		statefile << allps[temps[i]] << (i == (temps.size() - 1) ? "\n" : "\t");
	}
	statefile << "path for C\n";
	for (int i = 0; i < (int)temps.size(); i++) {
		statefile << allpc[temps[i]] << (i == (temps.size() - 1) ? "\n" : "\t");
	}
	statefile << "sigma " << allsigma[index] << '\n';
	statefile << "covariance matrix\n";
	for (int i = 0; i < (int)temps.size(); i++) {
		for (int j = 0; j <= i; j++) {
			statefile << allC[temps[i]][temps[j]] << (j == i ? "\n" : "\t");
		}
	}
	statefile.close();
}

CBOG_CBD::~CBOG_CBD() {
	delete[] gbest;
	for (int i = 0; i < DIM; i++) {
		delete[] allC[i];
	}
	delete[] allC;
	delete[] allmean;
	delete[] allps;
	delete[] allpc;
	for (int i = (int)optimizers.size() - 1; i >= 0; i--) {
		delete optimizers[i];
	}
	delete[] contribution;
	result.close();
	endtime = clock();
	ofstream timefile("CBCCOtimefile.txt", ios::app);
	timefile << fp->getID() << "," << (endtime - starttime) / CLOCKS_PER_SEC << endl;
	timefile.close();
}