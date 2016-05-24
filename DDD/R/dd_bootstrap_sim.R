dd_bootstrap_sim = function(brts, initparsopt = if(ddmodel < 5) {c(0.2,0.1,2*(length(brts) + missnumspec))} else {c(0.2,0.1,2*(length(brts) + missnumspec),0.01)}, idparsopt = 1:length(initparsopt), idparsfix = (1:(3 + (ddmodel == 5)))[-idparsopt], parsfix = (ddmodel < 5) * c(0.2,0.1,2*(length(brts) + missnumspec))[-idparsopt] + (ddmodel == 5) * c(0.2,0.1,2*(length(brts) + missnumspec),0)[-idparsopt], res = 10*(1+length(brts)+missnumspec), ddmodel = 1, missnumspec = 0, cond = 3, btorph = 1, soc = 2, tol = c(1E-3, 1E-4, 1E-6), maxiter = 1000 * round((1.25)^length(idparsopt)), endmc = 100, seed = 42)
                                               
{
# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# initparsopt contains initial parameter values
# - initparsopt[1] = b (= la_1 in ER2012) = speciation initiation rate
# - initparsopt[2] = mu_1 (= mu_g in ER2012) = extinction rate of good species
# - initparsopt[3] = la_1 (= la_2 in ER2012) = speciation completion rate
# - initparsopt[4] = mu_2 (= mu_i in ER2012) = extinction rate of incipient species
# res = resolution of the method; res should be larger than the total number of species
# missnumspec = number of missing species    
# cond = conditioning
# . cond = 0 conditioning on stem or clade age
# . cond = 1 conditioning on age and non-extinction of the phylogeny 
# . cond = 2 conditioning on age and on number of extant taxa
# soc = stem (1) or crown (2) age # Only works for crown age so far!
# btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# tol = tolerance in optimization
# - reltolx = relative tolerance of parameter values in optimization
# - reltolf = relative tolerance of function value in optimization
# - abstolx = absolute tolerance of parameter values in optimization
# maxiter = the maximum number of iterations in the optimization
# endmc = number of simulations in the bootstrap
# seed = seed for the simulations in the bootstrap

set.seed(seed)
cat('Finding the maximum likelihood estimates ...\n\n')
if(soc != 2)
{
   cat('This only works for clades starting from crown age!')
} else {
   empout = dd_ML(brts,initparsopt,idparsopt,idparsfix,parsfix,res,ddmodel,missnumspec,cond,btorph,soc,tol,maxiter)
   MLpars = as.numeric(c(empout[1:3]))
   empout = cbind(ntips = length(brts) + 1,empout)
   simout = NULL
   cat('Bootstrapping ...\n\n')
   trees = list()
   for(i in 1:endmc)
   {
      cat('Simulated data set',i,'out of',endmc,'\n')
      flush.console()
      simbrts = dd_sim(pars = c(MLpars), age = max(abs(brts)), ddmodel = ddmodel)
      trees[[i]] = brts2phylo(simbrts)
   }
   return(list(empout,trees))
}
}