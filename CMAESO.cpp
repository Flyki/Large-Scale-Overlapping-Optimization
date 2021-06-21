#include "CMAESO.h"

CMAESO::CMAESO(Benchmarks* fp, int id, vector<int> dims, int swarmsize, int seed) {
	this->fp = fp;
	this->id = id;
	this->novar = dims.size();
	stringstream ss;
	string idStr, theFuncStr;
	ss << id;
	ss >> idStr;
	ss.clear();
	ss << fp->getID();
	ss >> theFuncStr;
	parameterfilename = theFuncStr + "_" + idStr + "_cmaes_initials.par";
	restartfilename = theFuncStr + "_" + idStr + "_resume.par";
	writeParameters(novar);
	fitness = cmaes_init(&evo, novar, NULL, NULL, 0, swarmsize, parameterfilename.c_str());
	this->swarmsize = (int)cmaes_Get(&evo, "lambda");
	this->dims = dims;
	this->lowerBound = new double[2];
	this->lowerBound[0] = fp->getMinX();
	this->lowerBound[1] = fp->getMinX();
	this->upperBound = new double[2];
	this->upperBound[0] = fp->getMaxX();
	this->upperBound[1] = fp->getMaxX();
	cmaes_boundary_transformation_init(&boundaries, lowerBound, upperBound, 2);
	this->times = 2;
}

void CMAESO::writeParameters(int D) {
	ofstream parfile(parameterfilename.c_str());
	parfile << "N " << D << endl;
	parfile << "initialX 1:" << endl;
	parfile << " " << ((double)fp->getMaxX() + fp->getMinX()) / 2.0 << endl;
	parfile << "initialStandardDeviations 1:" << endl;
	parfile << " " << ((double)fp->getMaxX() - fp->getMinX()) / 3.0 << endl;
	parfile << "stopMaxFunEvals 1e299 " << endl;
	parfile << "stopMaxIter 1e299 " << endl;
	parfile << "stopTolFun 1e-100" << endl;
	parfile << "seed 0" << endl;
	parfile.close();
}

void CMAESO::writeParameters(int D, double* xmean) {
	ofstream parfile(parameterfilename.c_str());
	parfile << "N " << D << endl;
	parfile << "initialX " << D << ":" << endl;
	for (int i = 0; i < D; i++) {
		parfile << " " << xmean[i];
	}
	parfile << endl;
	parfile << "initialStandardDeviations 1:" << endl;
	parfile << " " << ((double)fp->getMaxX() - fp->getMinX()) / (3.0 * times) << endl;
	parfile << "stopMaxFunEvals 1e299 " << endl;
	parfile << "stopMaxIter 1e299 " << endl;
	parfile << "stopTolFun 1e-100" << endl;
	parfile << "seed 0" << endl;
	parfile.close();
}

void CMAESO::resume(vector<int> dims, int swarmsize) {
	cmaes_exit(&evo);
	this->dims = dims;
	this->novar = dims.size();
	writeParameters(dims.size());
	fitness = cmaes_init(&evo, dims.size(), NULL, NULL, 0, swarmsize, parameterfilename.c_str());
	this->swarmsize = (int)cmaes_Get(&evo, "lambda");
	char* resName = new char[restartfilename.length() + 1];
	strcpy_s(resName, restartfilename.length() + 1, restartfilename.c_str());
	cmaes_resume_distribution(&evo, resName);
	delete[] resName;
}

void CMAESO::restart(double* gbest) {
	double* xm = new double[novar];
	for (int i = 0; i < novar; i++) {
		xm[i] = gbest[dims[i]];
	}
	cmaes_exit(&evo);
	writeParameters(novar, xm);
	fitness = cmaes_init(&evo, novar, xm, NULL, 0, swarmsize, parameterfilename.c_str());
	delete[] xm;
	times++;
}

int CMAESO::calculateFitness(double* contextVector, int allDim, long int& fes) {
	double* temp = new double[allDim];
	memcpy(temp, contextVector, sizeof(double) * allDim);
	double* x_in_bounds = cmaes_NewDouble(novar);
	for (int i = 0; i < swarmsize; i++) {
		cmaes_boundary_transformation(&boundaries, swarm[i], x_in_bounds, novar);
		for (int j = 0; j < novar; j++) {
			temp[dims[j]] = x_in_bounds[j];
		}
		fitness[i] = fp->compute(temp);
		fes++;
	}
	curBest = 0;
	for (int i = 0; i < swarmsize; i++) {
		if (fitness[i] < fitness[curBest]) {
			curBest = i;
		}
	}
	free(x_in_bounds);
	delete[] temp;
	return curBest;
}

bool CMAESO::isFeasible(double* xx) {
	bool flag = true;
	for (int i = 0; i < novar; i++) {
		if (xx[i] > fp->getMaxX() || xx[i] < fp->getMinX()) {
			flag = false;
			break;
		}
	}
	return flag;
}

void CMAESO::generateSubsolution() {
	double* x_in_bounds = cmaes_NewDouble(novar);
	swarm = cmaes_SamplePopulation(&evo);
	for (int i = 0; i < swarmsize; i++) {
		cmaes_boundary_transformation(&boundaries, swarm[i], x_in_bounds, novar);
		while (isFeasible(x_in_bounds) == false) {
			swarm = cmaes_ReSampleSingle(&evo, i);
			cmaes_boundary_transformation(&boundaries, swarm[i], x_in_bounds, novar);
		}
	}
	free(x_in_bounds);
}

void CMAESO::updateTheDistribution() {
	cmaes_UpdateDistribution(&evo, fitness);
}

void CMAESO::updateContextVector(double* gbest, double& gbestf, int DIM) {
	if (fitness[curBest] < gbestf) {
		double* x_in_bounds = cmaes_NewDouble(novar);
		cmaes_boundary_transformation(&boundaries, swarm[curBest], x_in_bounds, novar);
		for (int i = 0; i < novar; i++) {
			gbest[dims[i]] = x_in_bounds[i];
		}
		gbestf = fitness[curBest];
		free(x_in_bounds);
	}
}

CMAESO::~CMAESO() {
	cmaes_exit(&evo);
	cmaes_boundary_transformation_exit(&boundaries);
	delete[] lowerBound;
	delete[] upperBound;
}