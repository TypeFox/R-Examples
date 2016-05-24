pbd_bootstrap_sim = function(brts, initparsopt = c(0.2,0.1,1), idparsopt = 1:length(initparsopt), idparsfix = NULL, parsfix = NULL, exteq = (length(initparsopt) < 4), parsfunc = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}), missnumspec = 0, cond = 1, btorph = 0, soc = 2, plotltt = 1, methode = "lsoda", n_low = 0, n_up = 0, tol = c(1E-4, 1E-4, 1E-6), maxiter = 1000 * round((1.25)^length(idparsopt)), endmc = 100, seed = 42)
{
# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# initparsopt contains initial parameter values
# - initparsopt[1] = b (= la_1 in ER2012) = speciation initiation rate
# - initparsopt[2] = mu_1 (= mu_g in ER2012) = extinction rate of good species
# - initparsopt[3] = la_1 (= la_2 in ER2012) = speciation completion rate
# - initparsopt[4] = mu_2 (= mu_i in ER2012) = extinction rate of incipient species
# exteq = incipient species have the same (1) or different (0) extinction rate as good species
# parsfunc = functions of parameters
# res = resolution of the method; res should be larger than the total number of species
# missnumspec = number of missing species    
# cond = conditioning
# . cond = 0 conditioning on stem or clade age
# . cond = 1 conditioning on age and non-extinction of the phylogeny 
# . cond = 2 conditioning on age and on number of extant taxa
# soc = stem (1) or crown (2) age
# btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# methode = method of the numerical integration; see package deSolve for details
# tol = tolerance in optimization
# - reltolx = relative tolerance of parameter values in optimization
# - reltolf = relative tolerance of function value in optimization
# - abstolx = absolute tolerance of parameter values in optimization
# maxiter = the maximum number of iterations in the optimization
# endmc = number of simulations in the bootstrap
# seed = seed for the simulations in the bootstrap

set.seed(seed)
cat('Finding the maximum likelihood estimates ...\n\n')
if(plotltt == 1)
{
    plot(c(-brts,0),c(soc:(length(brts)+soc-1),length(brts)+soc-1),log = 'y',type = 's',xlab = 'Time',ylab = 'Number of lineages')
}
empout = pbd_ML(brts,initparsopt,idparsopt,idparsfix,parsfix,exteq,parsfunc,missnumspec,cond,btorph,soc,methode,n_low,n_up,tol,maxiter)
MLpars = as.numeric(c(empout[1:4]))
exp_durspec = pbd_durspec_mean(c(MLpars[1],MLpars[3],MLpars[4]))
empout = cbind(ntips = length(brts) + 1,empout,exp_durspec)
simout = NULL
cat('Bootstrapping ...\n\n')
trees = list()
for(i in 1:endmc)
{
   cat('Simulated data set',i,'out of',endmc,'\n')
   flush.console()
   simbrts = pbd_sim_cpp(pars = MLpars, age = max(abs(brts)), soc = soc, plotltt = plotltt, methode = methode)
   trees[[i]] = brts2phylo(simbrts)
}
return(list(empout,trees))
}