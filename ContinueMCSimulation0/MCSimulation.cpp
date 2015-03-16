//============================================================================
// Name        : Monte Carlo Simulation
// Author      : amifan
// Version     : V1.0
//============================================================================

#include "Calculator.h"
#include "LogNormalDistribution.h"

using namespace std;

int main() {
	
	double tolerant = 0.11;
	double interval = 1E-10;
	double simulationTime = 0.00000005;
	unsigned int maxNumberOfParticles = 96000;

	ifstream wasser;
	ifstream glycerin;
	ifstream stickstoff;
	
	wasser.open ("0wasserDampfMasse");
	glycerin.open ("0glycerinDampfMasse");
	stickstoff.open ("0stickStoffDampfMasse");

	double wasserDampfMasse;
	double glycerinDampfMasse;
	double stickStoffDampfMasse;

	wasser >> wasserDampfMasse;
	glycerin >> glycerinDampfMasse;
	stickstoff >> stickStoffDampfMasse;

	cout<<"water: "<<wasser<<endl;
	cout<<"glycerin: "<<glycerin<<endl;
	cout<<"stickstoff: "<<stickstoff<<endl;
				
	
	/*wasserDampfMasse = 0.475694227;
	glycerinDampfMasse = 0.811322932;
	stickStoffDampfMasse = 0.14799376;*/

	// The average size of the particle in the initial state
	double x50 = 4E-9;
	// Sigma
	double sigma = 1.1;
	// Define the total number of all particles
	double totalNumberOfParticles = 16000;
	// Define the number of one kind of particle
	double maxNumberOfParticlesInOneClass = 100;
	// Volume size96000

	//double volumen = 1E-17;
	//cout<<"TEST1"<<endl;
	LogNormalDistribution logNormal = LogNormalDistribution(sigma, x50, tolerant, totalNumberOfParticles, maxNumberOfParticlesInOneClass, wasserDampfMasse, glycerinDampfMasse, stickStoffDampfMasse);
	//cout<<"TEST2"<<endl;

	// The following is to run simulation only with wachstum.
	// Calculator cal(logNormal, simulationTime, tolerant, maxNumberOfParticles, interval, wasserDampfMasse, glycerinDampfMasse, stickStoffDampfMasse);
	// cal.DoOnlyWachstum();


	// The following is to run simulation with all the processes.
	//Calculator cal(simulationTime, tolerant, maxNumberOfParticles, interval, wasserDampfMasse, glycerinDampfMasse, stickStoffDampfMasse);
	//cal.DoSimulation();
	cout<<"TEST3"<<endl;
	// The following is to run simulation with all the processes and with custom initialization.
	Calculator cal(logNormal, simulationTime, tolerant, maxNumberOfParticles, interval, wasserDampfMasse, glycerinDampfMasse, stickStoffDampfMasse);
	cout<<"PRE DO SIM"<<endl;
	cal.DoSimulation();
	cout<<"TEST4"<<endl;
	cout << "Finish!!" << endl;
}
