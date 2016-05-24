## Program Dingli Model Diagnosis. Using a mathematical model described by Dingli and coworkers, analyzes and produces probability distributions for the latency time to diagnosis from one leukemic stem cell.
## Most of the parameters and procedures are outlined in the following papers:
## Dingli D, Traulsen A, Pacheco JM (2007) Compartmental Architecture and Dynamics of Hematopoiesis. PLoS ONE 2(4): e345. doi:10.1371/journal.pone.0000345
## Dingli D, Traulsen A, Pacheco JM (2008) Chronic Myeloid Leukemia: Origin, Development, Response to Therapy, and Relapse. Clinical Leukemia, Vol. 2, No. 2, 133-139. 
## Lenaerts T, Pacheco JM, Traulsen A, and Dingli D. Tyrosine kinase inhibitor therapy can cure chronic myeloid leukemia without hitting leukemic stem cells. Haematologica 2010;95:900-907. doi:10.3324/haematol.2009.015271

## Copyright (C) 2011 Julian Landaw and Rainer Sachs
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; version 2.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License http://www.r-project.org/ for more details.
## Contact Julian Landaw, jlandaw@berkeley.edu or Rainer Sachs (sachs@math.berkeley.edu) for questions or comments.

## The main purpose of the following code is for running the function 'diagnosisSim', which finds  a random sample of latency times between the first leukemic stem cell and diagnosis according to the model by Dingli and coworkers). After running the script you can type "diagnosisSim()", which will output a list of the following:
## a vector of 10 latency times before diagnosis,
## the total number of simulations it required to get 10 diagnoses,
## and the number of diagnoses for which there still happened to be a leukemic stem cell.

## To change the number of diagnosis times to, e.g., 100 (instead of 10), simply type in "diagnosisSim(diagnoses = 100)". 


lam = 14^(1/4); 					## N_i/N_i-1, ratio of sizes between 2 adjacent compartments. approximates to ~1.93. the number of stem cells, N_0, is about 400
rrat = 1825^(1/32); 				## r_i/r_i-1, ratio of rates of differentiation between 2 adjacent compartments. approximates to ~1.26
eps = rrat*lam/(2*(rrat*lam-1)); 	## Probability of differentiation for healthy cells
epsCml = 0.72; 					## ... for leukemic cells
rStoc = rrat^(1:6); 				## Rates of replication for the stem cells and next 6 compartments (the stem cell and stochastic compartments)
rDeterm = rrat^(7:30); 				## Rates of replication for the next 25 compartments (the deterministic compartments)
stage1Cells = 1/((2*eps-1)*rrat/400); ## The number of healthy cells in the first stochastic compartment at equilibrium (number of cells coming in = number of cells going out) 

moran = function(regIni=399, cmlIni=1, period, timestep = 0.0025) 
## This function simulates the progression of the population of healthy/leukemic stem cells and determines how many healthy/leukemic cells will be entering the next compartment (due to differentiation)
## 'regIni' is the initial number of healthy stem cells
## 'cmlIni' initial number of leukemic stem cells 
## 'timestep' is in years 
## the average amount of time for a replication/export is 1/400 of a year. 

## 'regIni' and 'cmlIni' may be changed, but the default is 399 healthy stem cells and 1 leukemic stem cell. 
## The default 'timestep' is somewhat arbitrarily chosen to be 1/400 ? 

## This is a subroutine that is used in the main function 
{
	numSteps = floor(period/timestep); 			## how many timesteps occur in the given period of time.
	regCells = c(regIni,rep(0,times=numSteps)); 	## initializes the progression of healthy cells - the initial number followed by N empty time slots.
	cmlCells = c(cmlIni,rep(0,times=numSteps)); 	## initializes the progression of leukemic cells
	regFeed = rep(0,times=numSteps); 			## initializes the amount of healthy cells moving to the next compartment
	cmlFeed = regFeed; 						## initializes the amount of leukemic cells moving to the next compartment
	prob = timestep; 							## the probability any cell replicates happens to be numerically equivalent to the timestep
	ran = rbinom(numSteps,400,prob); 			## binomial distribution to determine how many moran processes occur at each timestep
	for (i in 2:(numSteps+1)) 
	{
		regCells[i] = regCells[i-1];
		cmlCells[i] = cmlCells[i-1];
		regFeed[i-1] = 0;
		cmlFeed[i-1] = 0;
		if (ran[i-1] > 0) ##if ran[i-1] = 0, no moran processes will occur during the time step - population remains stagnant and no differentiation
		{
			for (j in 1:ran[i-1]) ##ran[i-1], determined by the binomial distribution, is the number of moran steps during the time step 
			{
				r = runif(1);
				if (r < cmlCells[i]/(cmlCells[i] + regCells[i]))
				{
					##leukemic cell is exported to next compartment
					cmlCells[i] = cmlCells[i] - 1;
					cmlFeed[i-1] = cmlFeed[i-1] + 1; 
				}
				else
				{
					##healthy cell is exported to next compartment
					regCells[i] = regCells[i] - 1;
					regFeed[i-1] = regFeed[i-1] + 1; 
				}
				r = runif(1);
				if (r < cmlCells[i]/(cmlCells[i] + regCells[i]))
				{
					##leukemic cell divides into 2 leukemic cells
					cmlCells[i] = cmlCells[i] + 1;
				}
				else
				{
					##healthy cell divides into 2 healthy cells
					regCells[i] = regCells[i] + 1;
				}
			}
		}
	}
	return(list(regCells,cmlCells,regFeed,cmlFeed,timestep)); 
}

stocs = function(regIni = round(stage1Cells*lam^(0:5)), cmlIni = rep(0,times=6), regFeed, cmlFeed, timestep = 0.0025)
## This function ('stocs' stands for 'stochastics'), with given parameters, simulates stochastically the progression of the next 6 compartments after the stem cell compartment, as well as determining how many cells will be entering the next compartment (due to differentiation).
## 'regIni' is a vector of the initial number of healthy cells in each compartment
## 'cmlIni' the initial number of leukemic cells
## 'timestep' is the length of time between readings, and is in years 
## 'regFeed' and 'cmlFeed' are the amounts of healthy and leukemic cells entering the first stochastic compartment at each time step (determined by the function moran).
{
	numSteps = length(regFeed);
	regFeedNext = rep(0,times=numSteps); 				## initializes what will be fed to the first deterministic compartment
	cmlFeedNext = regFeedNext;
	regCells = c(regIni,rep(0,times=6*numSteps));
	dim(regCells) = c(6,numSteps+1); 				## matrix where the first column vector is the initial number of healthy cells in each compartment, followed by N empty column vectors
	cmlCells = c(cmlIni,rep(0,times=6*numSteps));
	dim(cmlCells) = c(6,numSteps+1); 				## first column vector is initial number of leukemic cells in each compartment
	prob = timestep*rStoc; 							## vector of six probabilities where the ith probability is the probability a cell in the ith compartment will divide
	for (i in 2:(numSteps+1))
	{
		regCells[1:6,i] = regCells[1:6,(i-1)];
		cmlCells[1:6,i] = cmlCells[1:6,(i-1)];
		regCells[1,i] = regCells[1,i] + regFeed[i-1]; ##first compartment gets fed the number of exported cells from the moran processes
		cmlCells[1,i] = cmlCells[1,i] + cmlFeed[i-1];
		for (j in 1:5)
		{
			p = prob[j];
			regTotalReps = rbinom(1,regCells[j,i-1],p); ## binomial distribution with p the probability and regCells[j,i-1] the number of 'trials'
			## determines how many cells divide
			cmlTotalReps = rbinom(1,cmlCells[j,i-1],p);
			regDifferentiates = rbinom(1,regTotalReps,eps); ## given how many divisions there are, binomial distribution determines how many differentiate 
			
			cmldifferentiates = rbinom(1,cmlTotalReps,epsCml);
			regReplicates = regTotalReps - regDifferentiates;
			cmlReplicates = cmlTotalReps - cmldifferentiates;
			regCells[j,i] = regCells[j,i] + regReplicates - regDifferentiates;
			cmlCells[j,i] = cmlCells[j,i] + cmlReplicates - cmldifferentiates;
			regCells[j+1,i] = regCells[j+1,i] + 2*regDifferentiates;
			cmlCells[j+1,i] = cmlCells[j+1,i] + 2*cmldifferentiates;
		}
		##6th compartment is treated the same, only here the cells are exported to the 7th compartment (via regFeedNext and cmlFeedNext), which is governed by the function 'determs'
		regTotalReps = rbinom(1,regCells[6,i-1],prob[6]);
		cmlTotalReps = rbinom(1,cmlCells[6,i-1],prob[6]);
		regDifferentiates = rbinom(1,regTotalReps,eps);
		cmldifferentiates = rbinom(1,cmlTotalReps,epsCml);
		regReplicates = regTotalReps - regDifferentiates;
		cmlReplicates = cmlTotalReps - cmldifferentiates;
		regCells[6,i] = regCells[6,i] + regReplicates - regDifferentiates;
		cmlCells[6,i] = cmlCells[6,i] + cmlReplicates - cmldifferentiates;
		regFeedNext[i-1] = 2*regDifferentiates;
		cmlFeedNext[i-1] = 2*cmldifferentiates; 
	} 
	return(list(regCells,cmlCells,regFeedNext,cmlFeedNext,timestep));
}

diffEq = function(y,params)
## A helper function, in which the 24 deterministic compartments are governed by the differential equation y' = diffEq(y, params)
## 'params' is a list of parameters supplying the probability of differentiation (eps or epsCml) and the 24 rates of replication (rDeterm)
## The sole purpose of this function is for solving the system of differential equations for the deterministic compartments (used in the next function, 'determs')
{
	primes = rep(0,times=24);
	e = params[[1]];
	r = params[[2]];
	primes[1] = -(2*e-1)*r[1]*y[1];
	primes[2:24] = 2*e*r[1:23]*y[1:23] - (2*e-1)*r[2:24]*y[2:24];
	return(primes);
}

determs = function(regIni = c(round(stage1Cells*lam^(6:29)),0), cmlIni = rep(0,times=25), regFeed, cmlFeed, timestep = 0.0025)
## This function ('determs' stands for 'deterministic'), with given parameters, simulates deterministically the progression of the last 25 compartments after the 'stochastic' compartments (governed by the function 'stocs'), and simulates how many mature cells will be produced at each timestep. 
## 'regIni' is a vector of the initial number of healthy cells in each compartment 
## 'cmlIni' the initial number of leukemic cells 
## 'timestep' is the length of time between readings, and is in years 
## 'regFeed' and 'cmlFeed' are the amounts of healthy and leukemic cells entering the first of the compartments at each time step (determined by the function 'stocs').
{
	numSteps = length(regFeed);
	regParams = list(eps,rDeterm);
	cmlParams = list(epsCml,rDeterm); ## Parameters used for the function 'diffEq' for the system of differential equations
	regCells = c(regIni,rep(0,times=25*numSteps));
	dim(regCells) = c(25,numSteps+1);
	cmlCells = c(cmlIni,rep(0,times=25*numSteps));
	dim(cmlCells) = c(25,numSteps+1);
	for (i in 2:(numSteps+1))
	{
		k1reg = timestep*diffEq(regCells[1:24,i-1],regParams); ## Runge-Kutta method for solving the system for the healthy cells
		k2reg = timestep*diffEq(regCells[1:24,i-1]+k1reg/2,regParams);
		k3reg = timestep*diffEq(regCells[1:24,i-1]+k2reg/2,regParams);
		k4reg = timestep*diffEq(regCells[1:24,i-1]+k3reg,regParams);
		regCells[1:24,i] = regCells[1:24,i-1] + 1/6*(k1reg + 2*k2reg + 2*k3reg + k4reg);
		regCells[1,i] = regCells[1,i] + regFeed[i-1];
		
		k1cml = timestep*diffEq(cmlCells[1:24,i-1],cmlParams); ## RK method for the leukemic cells
		k2cml = timestep*diffEq(cmlCells[1:24,i-1]+k1cml/2,cmlParams);
		k3cml = timestep*diffEq(cmlCells[1:24,i-1]+k2cml/2,cmlParams);
		k4cml = timestep*diffEq(cmlCells[1:24,i-1]+k3cml,cmlParams);
		cmlCells[1:24,i] = cmlCells[1:24,i-1] + 1/6*(k1cml + 2*k2cml + 2*k3cml + k4cml);
		cmlCells[1,i] = cmlCells[1,i] + cmlFeed[i-1];
	}
	regCells[25,] = timestep*(2*eps)*rDeterm[24]*regCells[24,];
	cmlCells[25,] = timestep*(2*epsCml)*rDeterm[24]*cmlCells[24,];  
	return(list(regCells,cmlCells,timestep));
}

simulate = function(regIni = c(399,round(stage1Cells*lam^(0:29)),0),cmlIni = c(1,rep(0,times=31)),period,timestep = 0.0025)
## A very simple function that runs the functions 'moran', 'stocs', and 'determs' to follow the progression of all the compartments together
## 'regIni' is a vector of the initial number of healthy cells in each compartment 
## 'cmlIni' the initial number of leukemic cells 
## 'period' is the length of time for which to simulate an individual
## 'timestep' is the length of time between readings, and is in years
{
	regStemCells = regIni[1];
	regStocCells = regIni[2:7];
	regDetermCells = regIni[8:32];
	cmlStemCells = cmlIni[1];
	cmlStocCells = cmlIni[2:7];
	cmlDetermCells = cmlIni[8:32];
	p = period;
	t = timestep;
	moranCells = moran(regIni = regStemCells,cmlIni = cmlStemCells,period = p,timestep = t);
	stocCells = stocs(regIni = regStocCells,cmlIni = cmlStocCells,regFeed = moranCells[[3]],cmlFeed = moranCells[[4]],timestep = t);
	determCells = determs(regIni = regDetermCells,cmlIni = cmlDetermCells,regFeed = stocCells[[3]],cmlFeed = stocCells[[4]],timestep = t);
	return(list(rbind(moranCells[[1]],stocCells[[1]],determCells[[1]]),rbind(moranCells[[2]],stocCells[[2]],determCells[[2]])));
}

diagnosisSim = function(regIni = c(399,round(stage1Cells*lam^(0:30)),0),cmlIni = c(1,rep(0,times=32)),period=50,timestep = 0.0025,diagnoses = 1000,stop=0)
## This function forms various simulations and outputs a list of the following: 
## a vector of various times for which CML diagnosis first comes 
## the total number of simulations it took to get so many latency times and 
## the number of those diagnoses in which there were still leukemic stem cells.

## 'regIni' is a vector of the initial number of healthy cells in each compartment. The default is the amount of cells of a healthy individual at equilibrium, minus one stem cell. 
## 'cmlIni' the initial number of leukemic cells. The default is one leukemic stem cell, and no leukemic cells in any other compartment.
## 'period' is the maximum length of time in years for which to simulate an individual. The default is 50, which is an extremely conservative maximum amount of time
## 'timestep' is the length of time between readings, and is in years
## 'diagnoses' is the number of diagnoses for which the user wants to output the latency times. The default is 10, and running the function will take around a minute of computer time.
## 'stop' allows the user to stop running the function once the number of simulations equals 'stop.' If stop = 0, as is the default, the function will run indefinitely until the number of diagnoses is met.

## If during a simulation it is clear an individual will not develop CML diagnosis, the simulation terminates and a new simulation begins
{
	diagnosisCount = 0;
	numTotalSims = 0;
	leukStemCellCount = 0;
	diagnosisTimes = rep(0,times = diagnoses);
	regular = regIni;
	leukemic = cmlIni; 
	while (diagnosisCount < diagnoses & (numTotalSims < stop | stop == 0))
	{
		numTotalSims = numTotalSims + 1;
		reg1 = regular;
		leuk1 = leukemic;
		for(i in 1:floor(period/timestep))
		{
			s = simulate(regIni = reg1,cmlIni = leuk1,period = timestep,timestep = timestep); ## simulates only the number of cells in the next time step
			regCells = s[[1]];
			cmlCells = s[[2]];
			r1 = regCells[,2];
			l1 = cmlCells[,2];
			if (leuk1[32] >= 10^12*365*timestep) ## diagnosis occurs - there are (on average) 10^12 mature leukemic cells produced per day.
			{
				diagnosisCount = diagnosisCount + 1;
				diagnosisTimes[diagnosisCount] = i*timestep;
				if (leuk1[1] > 0)
				{
					leukStemCellCount = leukStemCellCount + 1;
				}
				break;
			}
			else if (sum(leuk1[1:7]) == 0 & l1[32] < leuk1[32])   
			{
				## There are no leukemic cells in the stochastic compartments, and the number of mature leukemic cells is decreasing. 
				## This implies the mature leukemic cells will continue decreasing, thus no diagnosis can occur. To save time, the simulation terminates here.
				break;
			}
			reg1 = r1;
			leuk1 = l1;
		}
	}
	return(list(diagnosisTimes,numTotalSims,leukStemCellCount));
}
