dd_MS_ML = function(brtsM, brtsS, tsplit, initparsopt = c(0.5, 0.1, 2 * (1 + length(brtsM) + length(brtsS) + sum(missnumspec)),(tsplit + max(brtsS))/2), parsfix = NULL, idparsopt = c(1:3,6), idparsfix = NULL, idparsnoshift = (1:6)[c(-idparsopt,(-1)^(length(idparsfix) != 0) * idparsfix)], res = 10*(1 + length(c(brtsM,brtsS)) + sum(missnumspec)), ddmodel = 1.3, missnumspec = 0, cond = 0, soc = 2, tol = c(1E-3, 1E-4, 1E-6), maxiter = 1000 * round((1.25)^length(idparsopt)), changeloglikifnoconv = FALSE, optimmethod = 'subplex', methode = 'analytical')
{
# brtsM, brtsS = branching times of main clade and subclade (positive, from present to past)
# - max(brtsM) = crown age
# - min(brtsM,brtsS) = most recent branching time
# - tsplit = the branching time where the subclade branches off from the main clade
# - idparsopt contains the ids of the parameters to be optimized, e.g. to optimize la, mu, K, K2 and tshift idparsopt = c(1,2,3,6,7)
# - initparsopt contains the starting values of the parameters to be optimized
# - idparsfix contains the ids of the parameters that are fixed and must not be optimized
# - parsfix contains the values of the fixed parameters
# - idparsnoshift contains the ids of the parameters la2, mu2 and K2 that do not shift, i.e. that need to be set equal to la, mu and K
# - pars[1] = la_M = (initial) speciation rate in main clade
# - pars[2] = mu_M = extinction rate in main clade
# - pars[3] = K = carrying capacity
# - pars[4] = la_S = (initial) speciation rate in subclade
# - pars[5] = mu_S = extinction rate in subclade
# - pars[6] = t_d = time of decoupling
# - res = resolution of the method; res should be larger than the total number of species
# - ddmodel = diversity-dependent model,mode of diversity-dependence
#  . ddmodel == 1.3: linear dependence in speciation rate with parameter K'
#  . ddmodel == 2.3: exponential dependence in speciation rate with parameter x
# - missnumspec = number of missing species    
# - cond = conditioning:
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on non-extinction of the phylogeny NYI
# - tol = tolerance in optimization
#  . reltolx = relative tolerance of parameter values in optimization
#  . reltolf = relative tolerance of function value in optimization
#  . abstolx = absolute tolerance of parameter values in optimization
# - maxiter = the maximum number of iterations in the optimization
# - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
# - optimmethod = 'subplex' (current default) or 'simplex' (default of previous versions)
# - methode = the method used in the numerical solving of the set of the ode's

options(warn = -1)
brtsM = sort(abs(as.numeric(brtsM)),decreasing = TRUE)
brtsS = sort(abs(as.numeric(brtsS)),decreasing = TRUE)
if(cond == 1 & soc == 1)
{
   cat("Conditioning on survival of a clade with stem age currently not implemented.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K = -1, lambda_S = -1, mu_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
if(is.numeric(brtsM) == FALSE || is.numeric(brtsS) == FALSE)
{ 
   cat("The branching times should be numeric.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K = -1, lambda_S = -1, mu_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
idparsnoshift = sort(idparsnoshift)
idpars = sort(c(idparsopt,idparsfix,idparsnoshift))
if((sum(idpars == (1:6)) != 6) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
{
   cat("The parameters to be optimized, fixed and not shifted are incoherent.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K = -1, lambda_S = -1, mu_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
namepars = c("la_M","mu_M","K","la_S","mu_S","t_d")
if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
cat("You are optimizing",optstr,"\n")
if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
cat("You are fixing",fixstr,"\n")
if(length(namepars[idparsnoshift]) == 0) { noshiftstr = "anything" } else { noshiftstr = namepars[idparsnoshift] }
cat("You are not shifting",noshiftstr,"\n")
cat("Optimizing the likelihood - this may take a while.","\n")
flush.console()
trparsopt = initparsopt/(1 + initparsopt)
trparsopt[which(initparsopt == Inf)] = 1
trparsfix = parsfix/(1 + parsfix)
trparsfix[which(parsfix == Inf)] = 1
pars2 = c(res,ddmodel,cond,tsplit,0,soc,tol,maxiter)
optimpars = c(tol,maxiter)
initloglik = dd_MS_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,pars2 = pars2,brtsM = brtsM,brtsS = brtsS,missnumspec = missnumspec, methode = methode)
cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
flush.console()
if(initloglik == -Inf)
{
   cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K = -1, lambda_S = -1, mu_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
out = optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = dd_MS_loglik_choosepar,trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,pars2 = pars2,brtsM = brtsM,brtsS = brtsS,missnumspec = missnumspec,methode = methode)
if(out$conv > 0)
{
   cat("Optimization has not converged. Try again with different initial values.\n")
   out2 = data.frame(row.names = "results",lambda_1 = -1, mu_1 = -1, K = -1, lambda_2 = -1, mu_2 = -1, t_d = -1, loglik = -1, df = -1, conv = unlist(out$conv))
} else {
MLtrpars = as.numeric(unlist(out$par))
MLpars = MLtrpars/(1-MLtrpars)
MLpars1 = rep(0,6)
MLpars1[idparsopt] = MLpars
if(length(idparsfix) != 0) {MLpars1[idparsfix] = parsfix }
if(length(idparsnoshift) != 0) { MLpars1[idparsnoshift] = MLpars1[idparsnoshift - 3] }
if(MLpars1[3] > 10^7){MLpars1[3] = Inf}
ML = as.numeric(unlist(out$fvalues))
out2 = data.frame(row.names = "results",lambda_M = MLpars1[1],mu_M = MLpars1[2],K = MLpars1[3], lambda_S = MLpars1[4], mu_S = MLpars1[5], t_d = MLpars1[6], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
if(out2$conv != 0 & changeloglikifnoconv == T) { out2$loglik = -Inf }
s1 = sprintf('Maximum likelihood parameter estimates: %f %f %f %f %f %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6])
s2 = sprintf('Maximum loglikelihood: %f',ML)
cat("\n",s1,"\n",s2,"\n")
out$par = list(MLpars1)
}
}
}
}
}
invisible(out2)
}