#MCMC
BayesmConstant.keep = 1             #keep every keepth draw for MCMC routines
BayesmConstant.nprint = 100         #print the remaining time on every nprint'th draw
BayesmConstant.RRScaling = 2.38     #Roberts and Rosenthal optimal scaling constant
BayesmConstant.w = .1               #fractional likelihood weighting parameter

#Priors
BayesmConstant.A = .01              #scaling factor for the prior precision matrix
BayesmConstant.nuInc = 3            #Increment for nu
BayesmConstant.a = 5                #Dirichlet parameter for mixture models
BayesmConstant.nu.e = 3             #degrees of freedom parameter for regression error variance prior
BayesmConstant.nu = 3               #degrees of freedom parameter for Inverted Wishart prior
BayesmConstant.agammaprior = .5     #Gamma prior parameter
BayesmConstant.bgammaprior = .1     #Gamma prior parameter

#DP
BayesmConstant.DPalimdef=c(.01,10)  #defines support of 'a' distribution
BayesmConstant.DPnulimdef=c(.01,3)  #defines support of nu distribution
BayesmConstant.DPvlimdef=c(.1,4)    #defines support of v distribution
BayesmConstant.DPIstarmin = 1       #expected number of components at lower bound of support of alpha
BayesmConstant.DPpower = .8         #power parameter for alpha prior
BayesmConstant.DPalpha = 1.0        #intitalized value for alpha draws
BayesmConstant.DPmaxuniq = 200      #storage constraint on the number of unique components
BayesmConstant.DPSCALE = TRUE       #should data be scaled by mean,std deviation before posterior draws
BayesmConstant.DPgridsize = 20      #number of discrete points for hyperparameter priors

#Mathematical Constants
BayesmConstant.gamma = .5772156649015328606

#BayesBLP
BayesmConstant.BLPVOmega = matrix(c(1,0.5,0.5,1),2,2)  #IW prior parameter of correlated shocks in IV bayesBLP
BayesmConstant.BLPtol = 1e-6