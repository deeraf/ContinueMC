/*
 * LogNormalDistribution.cpp
 *
 *  Created on: Apr 6, 2013
 *      Author: amifan
 */

#include "LogNormalDistribution.h"
#include "Particle.h"
#include "ParticleClassVector.h"
#include "ParticleSystem.h"
#include "Logger.h"
#include "math.h"
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

LogNormalDistribution::LogNormalDistribution(double s, double x, double t, double N, double c,double wasserDampMasse,
		double glycerinDampMasse, double stickstoffDampMasse) {
	sigma = s;
	x50 = x;
	tolerant = t;
	totalNumberOfParticles = N;
	maxNumberOfParticlesInOneClass = c;
	ps = ParticleSystem(wasserDampMasse, glycerinDampMasse, stickstoffDampMasse);
}

double LogNormalDistribution::q(double x, double x50) {
	double result = 0.0;

	double a = sqrt(2 * M_PI) * log(sigma);
	double b = log(x / x50) / log(sigma);

	if (x != 0) {
		result = 1 / a * 1 / x * exp(-0.5 * pow(b, 2));
	} else {
		result = 0;
	}

	return result;
}
ParticleClassVector LogNormalDistribution::Initialize(const char *fname1, const char * fname2) {
	ParticleClassVector pcv;
	ifstream f1;
	ifstream f2;
	int num1=0, num2=0;
	
	f2.open("0classSize.txt");
	if(f2.fail()){
	      cout << "error2" << endl;
	}
	int n = std::distance(std::istream_iterator<double>(f2), std::istream_iterator<double>());
	int nopArray[n];
	double diamArray[n];
	f2.close();

	//cout<<"n="<<n<<endl;

	f1.open("0_MCOutput.txt");
	if(f1.fail()){
	      cout << "error" << endl;
	}
	while (!f1.eof() && num1<n){
		f1 >> nopArray[num1];
		//cout<<"nopArray: "<<nopArray[num1]<<endl;
		num1++;
	}
	f1.close();
	f2.open("0classSize.txt");
	if(f2.fail()){
	      cout << "error2" << endl;
	}
	
	while (!f2.eof() && num2<n){
		f2 >> diamArray[num2];
		//cout<<"diamArray: "<<diamArray[num2]<<endl;
		num2++;
	}
	f2.close();
/*	
    //Sort 
     int temp;
     for(int i2=0; i2<=n-1; i2++){
    	 for(int j=0; j<n-1; j++){  
           if(nopArray[j]>nopArray[j+1]){
        	   temp=nopArray[j];
        	   nopArray[j]=nopArray[j+1];
        	   nopArray[j+1]=temp;        
           }
    	 }         
     } 	
	cout<<"sorted"<<endl;    
	for(int i=0; i<n; i++){
		//f1>>nopArray[i];
		cout<<" "<<nopArray[i];
		//f2>>diamArray[i];
		//cout<<"1diamArray: "<<diamArray[i]<<endl;
	} */
/**/
	/*
	int temp, flag=1;
	for(int i = 1; (i <= n)&& flag; i++){
		flag = 0;
	          for (int j=0; j < (n -1); j++){
	               if (nopArray[j+1] < nopArray[j]){ 
	                    temp = nopArray[j];             // swap elements
	                    nopArray[j] = nopArray[j+1];
	                    nopArray[j+1] = temp;
	                    flag = 1;
	               }
	          }
	     }
						
	for (int i; i<n; i++){
		cout<<"1nopArray: "<<nopArray[i]<<endl;
	}
*/		
	for(int i=0; i < n; i++){
		pcv.InsertNewParticleClass();
		//cout<<"value of n: "<<n<< " \tvalue of i:"<<i<<endl;
		for (int j=0; j<nopArray[i]; j++){
		//	cout <<"j; "<<j<<"dimarray: "<<diamArray[i]<<endl;
			//cout<<"j: "<<j <<" nopArray: "<<nopArray[i]<<endl;
			pcv.AddParticle(i, diamArray[i]); //pcv.AddParticle(i, diamArray[i]);
			//cout<<"particle added with success"<<endl;
		}
		//cout<<"TEST12"<<endl;
		pcv.UpdateBoundary();
		//cout<<"TEST13"<<endl;
	}
	//cout<<"TEST10"<<endl;
	return pcv;	
}

ParticleClassVector LogNormalDistribution::Initialize() {
	ParticleClassVector pcv;


	double deltaX = 1E-15;
	double xi = 0;
	double interval = 1 / totalNumberOfParticles;
	double currentR = 0;
	double nextR = 0;
	double *xiArray = (double*) malloc(sizeof(double) * totalNumberOfParticles);
	double diff1;
	double diff2;


	for (int i = 0; i < totalNumberOfParticles; i++) {
		nextR = (i + 0.5) * interval;

		while (currentR < nextR) {
			xi += deltaX;
			currentR += q(xi - deltaX * 0.5, x50) * deltaX;
		}

		xiArray[i] = xi;
		cout << "deltaA: " << i << " xi " << xi << endl;
	}

	do {
		int countOfParticleInClass = 0;
		int index = 0;

		pcv = ParticleClassVector(maxNumberOfParticlesInOneClass, tolerant);

		pcv.InsertNewParticleClass();

		for (int j = 0; j < totalNumberOfParticles; j++) {
			if (countOfParticleInClass < maxNumberOfParticlesInOneClass) {
				pcv.AddParticle(index, xiArray[j]);
				countOfParticleInClass++;
			} else {
				pcv.InsertNewParticleClass();
				index++;
				double durchmesser = xiArray[j];
				unsigned int numberOfWasser = ps.CalWasserNumber(durchmesser);
				unsigned int numberOfGlycerin = ps.CalGlycerinNumber(durchmesser);
				double kritischMolWasser = ps.CalkritischeMolWasser();
				double kritischMolGlycerin = ps.CalkritischeMolGlyzerin();
				double zusammensetzung = ps.CalKritischGlycerinZusammensetzung();
//				double zusammensetzung = 0.0;

				pcv.AddParticle(index, xiArray[j], zusammensetzung, kritischMolWasser, kritischMolGlycerin, numberOfWasser, numberOfGlycerin);
//				pcv.AddParticle(index, xiArray[j]);
				countOfParticleInClass = 1;
			}
		}

		for (int k = 0; k < pcv.GetNumberOfParticleClasses(); k++) {
			pcv.GetParticleClassAt(k).CalAverageDiameter();
		}

		pcv.UpdateBoundary();

		diff1 = (pcv.GetParticleClassAt(0).GetMaxDiameter() - pcv.GetParticleClassAt(0).GetMinDiameter()) / pcv.GetParticleClassAt(0).GetAverageDiameter();
		diff2 = (pcv.GetParticleClassAt(pcv.GetNumberOfParticleClasses() - 1).GetMaxDiameter() - pcv.GetParticleClassAt(pcv.GetNumberOfParticleClasses() - 1).GetMinDiameter())
				/ pcv.GetParticleClassAt(pcv.GetNumberOfParticleClasses() - 1).GetAverageDiameter();

		maxNumberOfParticlesInOneClass -= 1;
	} while (diff1 > tolerant || diff2 > tolerant);

	for (int l = 0; l < pcv.GetNumberOfParticleClasses(); l++) {
		for (int m = 0; m < pcv.GetParticleClassAt(l).GetSize(); m++) {
			pcv.GetParticleClassAt(l).GetParticleAt(m).SetWasserMol(ps.CalkritischeMolWasser());
			pcv.GetParticleClassAt(l).GetParticleAt(m).SetGlycerinMol(ps.CalkritischeMolGlyzerin());
			cout<<pcv.GetParticleClassAt(l).GetParticleAt(m).GetWasserMol()<<"  "<<pcv.GetParticleClassAt(l).GetParticleAt(m).GetGlycerinMol()<<endl;
		}
	}

	//Logger log;
	//log.GetMaxNumberOfParticlesAccordingTotTolerant(++maxNumberOfParticlesInOneClass);

	return pcv;
}
