dd_KI_ML = function(brtsM, brtsS, tsplit, initparsopt = c(0.5,0.1,2*(1 + length(brtsM) + missnumspec[1]),2*(1 + length(brtsS) + missnumspec[length(missnumspec)]),(tsplit + max(brtsS))/2), parsfix = NULL, idparsopt = c(1:3,6:7), idparsfix = NULL, idparsnoshift = (1:7)[c(-idparsopt,(-1)^(length(idparsfix) != 0) * idparsfix)], res = 10*(1 + length(c(brtsM,brtsS)) + sum(missnumspec)), ddmodel = 1, missnumspec = 0, cond = 1, soc = 2, tol = c(1E-3, 1E-4, 1E-6), maxiter = 1000 * round((1.25)^length(idparsopt)), changeloglikifnoconv = FALSE, optimmethod = 'subplex',methode = 'analytical')
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
# - pars[3] = K_M = carrying capacity in main clade
# - pars[4] = la_S = (initial) speciation rate in subclade
# - pars[5] = mu_S = extinction rate in subclade
# - pars[6] = K_S = carrying capacity in subclade
# - pars[7] = t_d = time of decoupling
# - res = resolution of the method; res should be larger than the total number of species
# - ddmodel = diversity-dependent model,mode of diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate with parameter K
#  . ddmodel == 1.3: linear dependence in speciation rate with parameter K'
#  . ddmodel == 2 : exponential dependence in speciation rate
#  . ddmodel == 2.1: variant with offset at infinity
#  . ddmodel == 2.2: 1/n dependence in speciation rate
#  . ddmodel == 2.3: exponential dependence in speciation rate with parameter x
#  . ddmodel == 3 : linear dependence in extinction rate
#  . ddmodel == 4 : exponential dependence in extinction rate
#  . ddmodel == 4.1: variant with offset at infinity
#  . ddmodel == 4.2: 1/n dependence in speciation rate
# - missnumspec = number of missing species    
# - cond = conditioning:
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on non-extinction of the phylogeny
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
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K_M = -1, lambda_S = -1, mu_S = -1, K_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
if(is.numeric(brtsM) == FALSE || is.numeric(brtsS) == FALSE)
{ 
   cat("The branching times should be numeric.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K_M = -1, lambda_S = -1, mu_S = -1, K_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
idparsnoshift = sort(idparsnoshift)
idpars = sort(c(idparsopt,idparsfix,idparsnoshift))
if((sum(idpars == (1:7)) != 7) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
{
   cat("The parameters to be optimized, fixed and not shifted are incoherent.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K_M = -1, lambda_S = -1, mu_S = -1, K_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {
namepars = c("la_M","mu_M","K_M","la_S","mu_S","K_S","t_d")
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
initloglik = dd_KI_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,pars2 = pars2,brtsM = brtsM,brtsS = brtsS,missnumspec = missnumspec,methode = methode)
cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
flush.console()
if(initloglik == -Inf)
{
   cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
   out2 = data.frame(row.names = "results",lambda_M = -1, mu_M = -1, K_M = -1, lambda_S = -1, mu_S = -1, K_S = -1, t_d = -1, loglik = -1, df = -1, conv = -1)
} else {

#code up to DDD v1.6: out = optimx2(trparsopt,dd_KI_loglik_choosepar,hess=NULL,method = "Nelder-Mead",hessian = FALSE,control = list(maximize = TRUE,abstol = pars2[8],reltol = pars2[7],trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,brtsM = brtsM, brtsS = brtsS, pars2 = pars2, missnumspec = missnumspec)
out = optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = dd_KI_loglik_choosepar,trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,pars2 = pars2,brtsM = brtsM,brtsS = brtsS,missnumspec = missnumspec, methode = methode)
if(out$conv > 0)
{
   cat("Optimization has not converged. Try again with different initial values.\n")
   out2 = data.frame(row.names = "results",lambda_1 = -1, mu_1 = -1, K_1 = -1, lambda_2 = -1, mu_2 = -1, K_2 = -1, t_d = -1, loglik = -1, df = -1, conv = unlist(out$conv))
} else {
MLtrpars = as.numeric(unlist(out$par))
MLpars = MLtrpars/(1-MLtrpars)
MLpars1 = rep(0,7)
MLpars1[idparsopt] = MLpars
if(length(idparsfix) != 0) {MLpars1[idparsfix] = parsfix }
if(length(idparsnoshift) != 0) { MLpars1[idparsnoshift] = MLpars1[idparsnoshift - 3] }
if(MLpars1[3] > 10^7){MLpars1[3] = Inf}
if(MLpars1[6] > 10^7){MLpars1[6] = Inf}
ML = as.numeric(unlist(out$fvalues))
out2 = data.frame(row.names = "results",lambda_M = MLpars1[1],mu_M = MLpars1[2],K_M = MLpars1[3], lambda_S = MLpars1[4], mu_S = MLpars1[5], K_S = MLpars1[6], t_d = MLpars1[7], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
if(out2$conv != 0 & changeloglikifnoconv == T) { out2$loglik = -Inf }
s1 = sprintf('Maximum likelihood parameter estimates: %f %f %f %f %f %f %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7])
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