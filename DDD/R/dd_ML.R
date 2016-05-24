dd_ML = function(brts, initparsopt = if(ddmodel < 5) {c(0.2,0.1,2*(length(brts) + missnumspec))} else {c(0.2,0.1,2*(length(brts) + missnumspec),0.01)}, idparsopt = 1:length(initparsopt), idparsfix = (1:(3 + (ddmodel == 5)))[-idparsopt], parsfix = (ddmodel < 5) * c(0.2,0.1,2*(length(brts) + missnumspec))[-idparsopt] + (ddmodel == 5) * c(0.2,0.1,2*(length(brts) + missnumspec),0)[-idparsopt], res = 10*(1+length(brts)+missnumspec), ddmodel = 1, missnumspec = 0, cond = 1, btorph = 1, soc = 2, tol = c(1E-3, 1E-4, 1E-6), maxiter = 1000 * round((1.25)^length(idparsopt)), changeloglikifnoconv = FALSE, optimmethod = 'subplex', methode = 'analytical')
{
# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - initpars[1] = la = (initial) speciation rate
# - initpars[2] = mu = extinction rate
# - initpars[3] = K = carrying capacity
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
#  . ddmodel == 5 : linear dependence in speciation rate and in extinction rate
# - missnumspec = number of missing species    
# - cond = conditioning
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on non-extinction of the phylogeny
#  . cond == 2 : conditioning on non-extinction of the phylogeny and on the total number of extant taxa (including missing species)
#  . cond == 3 : conditioning on the total number of extant taxa (including missing species)
# - btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# - tol = tolerance in optimization
#  . reltolx = relative tolerance of parameter values in optimization
#  . reltolf = relative tolerance of function value in optimization
#  . abstolx = absolute tolerance of parameter values in optimization
# - maxiter = the maximum number of iterations in the optimization
# - changeloglikifnoconv = if T the loglik will be set to -Inf if ML does not converge
# - optimmethod = 'subplex' (current default) or 'simplex' (default of previous versions)
# - methode = the method used in the numerical solving of the set of the ode's

options(warn=-1)
brts = sort(abs(as.numeric(brts)),decreasing = TRUE)
if(is.numeric(brts) == FALSE)
{
   cat("The branching times should be numeric.\n")
   out2 = data.frame(lambda = -1,mu = -1,K = -1, loglik = -1, df = -1, conv = -1)
   if(ddmodel == 5) {out2 = data.frame(lambda = -1,mu = -1,K = -1, r = -1, loglik = -1, df = -1, conv = -1)}
} else {
idpars = sort(c(idparsopt,idparsfix))
if((sum(idpars == (1:3)) != 3) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
{
   cat("The parameters to be optimized and/or fixed are incoherent.\n")
   out2 = data.frame(lambda = -1,mu = -1,K = -1, loglik = -1, df = -1, conv = -1)
   if(ddmodel == 5) {out2 = data.frame(lambda = -1,mu = -1,K = -1, r = -1, loglik = -1, df = -1, conv = -1)}
} else {
namepars = c("lambda","mu","K")
if(ddmodel == 5) {namepars = namepars = c("lambda","mu","K","r")}
if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
cat("You are optimizing",optstr,"\n")
if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
cat("You are fixing",fixstr,"\n")
cat("Optimizing the likelihood - this may take a while.","\n")
flush.console()
trparsopt = initparsopt/(1 + initparsopt)
trparsopt[which(initparsopt == Inf)] = 1
trparsfix = parsfix/(1 + parsfix)
trparsfix[which(parsfix == Inf)] = 1
pars2 = c(res,ddmodel,cond,btorph,0,soc,tol,maxiter)
optimpars = c(tol,maxiter)
initloglik = dd_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,pars2 = pars2,brts = brts,missnumspec = missnumspec, methode = methode)
cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
flush.console()
if(initloglik == -Inf)
{
   cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
   out2 = data.frame(lambda = -1,mu = -1,K = -1, loglik = -1, df = -1, conv = -1)
   if(ddmodel == 5) {out2 = data.frame(lambda = -1,mu = -1,K = -1, r = -1, loglik = -1, df = -1, conv = -1)}
} else {
#code up to DDD v1.6: out = optimx2(trparsopt,dd_loglik_choosepar,hess=NULL,method = "Nelder-Mead",hessian = FALSE,control = list(maximize = TRUE,abstol = pars2[8],reltol = pars2[7],trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,brts = brts, pars2 = pars2,missnumspec = missnumspec)
#out = dd_simplex(trparsopt,idparsopt,trparsfix,idparsfix,pars2,brts,missnumspec)
out = optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = dd_loglik_choosepar,trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,pars2 = pars2,brts = brts, missnumspec = missnumspec, methode = methode)
if(out$conv != 0)
{
   cat("Optimization has not converged. Try again with different initial values.\n")
   out2 = data.frame(lambda = -1,mu = -1,K = -1, loglik = -1, df = -1, conv = unlist(out$conv))
   if(ddmodel == 5) {out2 = data.frame(lambda = -1,mu = -1,K = -1, r = -1, loglik = -1, df = -1, conv = unlist(out$conv))}
} else {
MLtrpars = as.numeric(unlist(out$par))
MLpars = MLtrpars/(1-MLtrpars)
MLpars1 = rep(0,3)
if(ddmodel == 5) {MLpars1 = rep(0,4)}
MLpars1[idparsopt] = MLpars
if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
if(MLpars1[3] > 10^7){MLpars1[3] = Inf}
ML = as.numeric(unlist(out$fvalues))
out2 = data.frame(lambda = MLpars1[1],mu = MLpars1[2],K = MLpars1[3], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
s1 = sprintf('Maximum likelihood parameter estimates: lambda: %f, mu: %f, K: %f',MLpars1[1],MLpars1[2],MLpars1[3])
if(ddmodel == 5)
{
   s1 = sprintf('%s, r: %f',s1,MLpars1[4])
   out2 = data.frame(lambda = MLpars1[1],mu = MLpars1[2],K = MLpars1[3], r = MLpars1[4], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))   
}
if(out2$conv != 0 & changeloglikifnoconv == T) { out2$loglik = -Inf }
s2 = sprintf('Maximum loglikelihood: %f',ML)
cat("\n",s1,"\n",s2,"\n")
}
}
}
}
invisible(out2)
}