rhierNegbinRw= function(Data, Prior, Mcmc) {
#   Revision History
#	  Sridhar Narayanan - 05/2005
#         P. Rossi 6/05
#         fixed error with nobs not specified and changed llnegbinFract 9/05
#         3/07 added classes
#         3/08 fixed fractional likelihood
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
#   Model
#       (y_i|lambda_i,alpha) ~ Negative Binomial(Mean = lambda_i, Overdispersion par = alpha)
#
#       ln(lambda_i) =  X_i * beta_i
#
#       beta_i = Delta'*z_i + nu_i
#               nu_i~N(0,Vbeta)
#
#   Priors
#       vec(Delta|Vbeta) ~ N(vec(Deltabar), Vbeta (x) (Adelta^-1))
#       Vbeta ~ Inv Wishart(nu, V)
#       alpha ~ Gamma(a,b) where mean = a/b and variance = a/(b^2)
#
#   Arguments
#       Data = list of regdata,Z 
#           regdata is a list of lists each list with members y, X
#              e.g. regdata[[i]]=list(y=y,X=X)
#              X has nvar columns including a first column of ones
#              Z is nreg=length(regdata) x nz with a first column of ones
#
#       Prior - list containing the prior parameters
#           Deltabar, Adelta - mean of Delta prior, inverse of variance covariance of Delta prior
#           nu, V - parameters of Vbeta prior
#           a, b - parameters of alpha prior
#
#       Mcmc - list containing
#           R is number of draws
#           keep is thinning parameter (def = 1)
#           nprint - print estimated time remaining on every nprint'th draw (def = 100)
#           s_beta - scaling parameter for beta RW (def = 2.93/sqrt(nvar))
#           s_alpha - scaling parameter for alpha RW (def = 2.93)
#           w - fractional weighting parameter (def = .1)
#           Vbeta0, Delta0 - initial guesses for parameters, if not supplied default values are used
#


#
# Definitions of functions used within rhierNegbinRw (but outside of Rcpp loop)
#
llnegbinR = function(par,X,y, nvar) {
# Computes the log-likelihood
    beta = par[1:nvar]
    alpha = exp(par[nvar+1])+1.0e-50
    mean=exp(X%*%beta)
    prob=alpha/(alpha+mean)
    prob=ifelse(prob<1.0e-100,1.0e-100,prob)
     out=dnbinom(y,size=alpha,prob=prob,log=TRUE)
     return(sum(out))
}

llnegbinFract = 
function(par,X,y,Xpooled, ypooled, w,wgt, nvar,lnalpha)  {
# Computes the fractional log-likelihood at the unit level
    theta = c(par,lnalpha)
    (1-w)*llnegbinR(theta,X,y,nvar) + w*wgt*llnegbinR(theta,Xpooled,ypooled, nvar) 
}

#
# Error Checking
#
if(missing(Data)) {pandterm("Requires Data argument -- list of regdata and (possibly) Z")}

if(is.null(Data$regdata)) {
    pandterm("Requires Data element regdata -- list of data for each unit : y and X")
}
regdata=Data$regdata
nreg = length(regdata)

if (is.null(Data$Z)) {
    cat("Z not specified - using a column of ones instead", fill = TRUE)
    Z = matrix(rep(1,nreg),ncol=1)
}
else {
    if (nrow(Data$Z) != nreg) {
        pandterm(paste("Nrow(Z) ", nrow(Z), "ne number units ",nreg))
    }
    else {
        Z = Data$Z
    }
}
nz = ncol(Z)

dimfun = function(l) {
    c(length(l$y),dim(l$X))
}
dims=sapply(regdata,dimfun)
dims = t(dims)
nvar = quantile(dims[,3],prob=0.5)
for (i in 1:nreg) {
        if (dims[i, 1] != dims[i, 2] || dims[i, 3] != nvar) {
            pandterm(paste("Bad Data dimensions for unit ", i, 
                " dims(y,X) =", dims[i, ]))
        }
}

ypooled = NULL
Xpooled = NULL
for (i in 1:nreg) {
    ypooled = c(ypooled,regdata[[i]]$y)
    Xpooled = rbind(Xpooled,regdata[[i]]$X)
}
nobs= length(ypooled)

nvar=ncol(Xpooled)
#
# check for prior elements
#
if(missing(Prior)) {
    Deltabar=matrix(rep(0,nvar*nz),nrow=nz) ; Adelta=BayesmConstant.A*diag(nz) ; nu=nvar+BayesmConstant.nuInc; V=nu*diag(nvar); a=0.5; b=0.1;
}
else {
    if(is.null(Prior$Deltabar)) {Deltabar=matrix(rep(0,nvar*nz),nrow=nz)} else {Deltabar=Prior$Deltabar}
    if(is.null(Prior$Adelta)) {Adelta=BayesmConstant.A*diag(nz)} else {Adelta=Prior$Adelta}
    if(is.null(Prior$nu)) {nu=nvar+BayesmConstant.nuInc} else {nu=Prior$nu}
    if(is.null(Prior$V)) {V=nu*diag(nvar)} else {V=Prior$V}
    if(is.null(Prior$a)) {a=BayesmConstant.agammaprior} else {a=Prior$a}
    if(is.null(Prior$b)) {b=BayesmConstant.bgammaprior} else {b=Prior$b}
}

if(sum(dim(Deltabar) == c(nz,nvar)) != 2) pandterm("Deltabar is of incorrect dimension")
if(sum(dim(Adelta)==c(nz,nz)) != 2) pandterm("Adelta is of incorrect dimension")
if(nu < nvar) pandterm("invalid nu value")
if(sum(dim(V)==c(nvar,nvar)) != 2) pandterm("V is of incorrect dimension")
if((length(a) != 1) | (a <=0)) pandterm("a should be a positive number")
if((length(b) != 1) | (b <=0)) pandterm("b should be a positive number")

#
# check for Mcmc 
#
if(missing(Mcmc)) pandterm("Requires Mcmc argument -- at least R")
if(is.null(Mcmc$R)) {pandterm("Requires element R of Mcmc")} else {R=Mcmc$R}
if(is.null(Mcmc$Vbeta0)) {Vbeta0=diag(nvar)} else {Vbeta0=Mcmc$Vbeta0}
if(sum(dim(Vbeta0) == c(nvar,nvar)) !=2) pandterm("Vbeta0 is not of dimension nvar")
if(is.null(Mcmc$Delta0)) {Delta0=matrix(rep(0,nz*nvar),nrow=nz)} else {Delta0=Mcmc$Delta0}
if(sum(dim(Delta0) == c(nz,nvar)) !=2) pandterm("Delta0 is not of dimension nvar by nz")
if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
  if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
if(is.null(Mcmc$s_alpha)) { s_alpha=BayesmConstant.RRScaling} 
    else {s_alpha= Mcmc$s_alpha }
if(is.null(Mcmc$s_beta)) { s_beta=BayesmConstant.RRScaling/sqrt(nvar)} 
    else {s_beta=Mcmc$s_beta }
if(is.null(Mcmc$w)) { w=BayesmConstant.w} 
    else {w = Mcmc$w}
# Wayne Taylor 12/2014 #############################################
if(is.null(Mcmc$alpha)) {fixalpha=FALSE} else {fixalpha=TRUE; alpha=Mcmc$alpha}
if(fixalpha & alpha<=0) pandterm("alpha is not positive")
###################################################################

# print out problem
#
cat(" ",fill=TRUE)
cat("Starting Random Walk Metropolis Sampler for Hierarchical Negative Binomial Regression",fill=TRUE)
cat("  ",nobs," obs; ",nvar," covariates (including the intercept); ",fill=TRUE)
cat("  ",nz," individual characteristics (including the intercept) ",fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parameters:",fill=TRUE)
cat("Deltabar",fill=TRUE)
print(Deltabar)
cat("Adelta",fill=TRUE)
print(Adelta)
cat("nu",fill=TRUE)
print(nu)
cat("V",fill=TRUE)
print(V)
cat("a",fill=TRUE)
print(a)
cat("b",fill=TRUE)
print(b)
cat(" ",fill=TRUE)
cat("MCMC Parameters:",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat("s_alpha = ",s_alpha,fill=TRUE)
cat("s_beta = ",s_beta,fill=TRUE)
cat("Fractional Likelihood Weight Parameter = ",w,fill=TRUE)
cat(" ",fill=TRUE)

par = rep(0,(nvar+1))
cat("initializing Metropolis candidate densities for ",nreg,"units ...",fill=TRUE)
fsh()
mle = optim(par,llnegbinR, X=Xpooled, y=ypooled, nvar=nvar, 
      method="L-BFGS-B", upper=c(Inf,Inf,Inf,log(100000000)), hessian=TRUE, control=list(fnscale=-1))
fsh()
beta_mle=mle$par[1:nvar]
alpha_mle = exp(mle$par[nvar+1])
varcovinv = -mle$hessian
Delta = Delta0
Beta = t(matrix(rep(beta_mle,nreg),ncol=nreg))
Vbetainv = chol2inv(chol(Vbeta0)) #Wayne: replaced "solve" function
Vbeta = Vbeta0
alpha = alpha_mle
alphacvar = s_alpha/varcovinv[nvar+1,nvar+1]
alphacroot = sqrt(alphacvar)
cat("beta_mle = ",beta_mle,fill=TRUE)
cat("alpha_mle = ",alpha_mle, fill = TRUE)
fsh()

hess_i=NULL
if(nobs > 1000){
  sind=sample(c(1:nobs),size=1000)
  ypooleds=ypooled[sind]
  Xpooleds=Xpooled[sind,]
  }
else{
	ypooleds=ypooled
	Xpooleds=Xpooled
}
# Find the individual candidate hessian
for (i in 1:nreg) {
    wgt = length(regdata[[i]]$y)/length(ypooleds)
    mle2 = optim(mle$par[1:nvar],llnegbinFract, X=regdata[[i]]$X, y=regdata[[i]]$y, Xpooled=Xpooleds, 
           ypooled=ypooleds, w=w,wgt=wgt, nvar=nvar, lnalpha=mle$par[nvar+1], 
           method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=0,reltol=1e-6))
    if (mle2$convergence==0)
        hess_i[[i]] = list(hess=-mle2$hessian)
    else
        hess_i[[i]] = diag(rep(1,nvar))
   if(i%%50 ==0) cat("  completed unit #",i,fill=TRUE)	
   fsh()
}

###################################################################
# Wayne Taylor
# 12/01/2014
###################################################################
if (fixalpha) {alpha=Mcmc$alpha}
rootA = chol(Vbetainv)
draws=rhierNegbinRw_rcpp_loop(regdata, hess_i, Z, Beta, Delta,
                              Deltabar, Adelta, nu, V, a, b,
                              R, keep, s_beta, alphacroot, 1000, rootA,
                              alpha, fixalpha)
###################################################################

attributes(draws$alphadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$alphadraw)$mcpar=c(1,R,keep)
attributes(draws$Deltadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$Deltadraw)$mcpar=c(1,R,keep)
attributes(draws$Vbetadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(draws$Vbetadraw)$mcpar=c(1,R,keep)
attributes(draws$Betadraw)$class=c("bayesm.hcoef")

return(draws)
}