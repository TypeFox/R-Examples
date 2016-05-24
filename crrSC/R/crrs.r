################################################################################
# Description: regression modeling of subdistribution hazards                  #
#              for stratified right censored data                              #
#              -- a modification to 'crr' in 'cmprsk' package                  #
# Reference:                                                                   #
#     -- R package: "cmprsk" by Robert Gray (directly used its subroutines)    #
#     -- manuscript: Competing Risks Regression for Stratified Data by         #
#                    Zhou, Latouche, Rocha, and Fine                           #
# Functions included:                                                          #
#     -- "crrs":  for fitting Ph subdistribution model for stratified data     #
#     -- "crrfsvs": to obtain logliklihood, score, and information             #
#     -- "crrfs": to obtain the logliklihood                                   #
#     -- "crrvvs": to obtain the variance estimators                           #
# Comments:                                                                    #
#     -- for highly stratified, bootstrap var estimate is recommended          #
# Author: Bingqing Zhou                                                        #
# Last updated: 05/04/2010  at 17:30                                           #
################################################################################

library(survival)

crrs = function(ftime, fstatus, cov1, cov2, strata, tf, failcode=1, cencode=0,
    ctype=1, subsets, na.action=na.omit, gtol=1e-6, maxiter=10,init) {
# arguments:
# ctype = 1 if estimating censoring dist within strata (regular stratification)
#       = 2 if estimating censoring dist across strata (highly stratification)
# ftime = vector of failure/censoring times
# fstatus = vector with a unique code for each failure type and a
#     separate code for censored observations
# cov1 = (nobs x ncovs) matrix of fixed covariates
# cov2 = matrix of covariates multiplied by functions of time;
#     if used, often these covariates would also appear in cov1,
#     to give a prop hazards effect plus a time interaction
# tf = functions of time.  A function that takes a vector of times as
#     an argument and returns a matrix whose jth column is the value of
#     the time function corresponding to the jth column of cov2 evaluated
#     at the input time vector.  At time tk, the
#     model includes the term cov2[,j]*tfs(tk)[,j] as a covariate.
# cengroup = vector with different values for each group with
#     a distinct censoring distribution (the censoring distribution
#     is estimated separately within these groups)
# failcode = code of fstatus that denotes the failure type of interest
# cencode = code of fstatus that denotes censored observations
# subset = logical vector length(ftime) indicating which cases to include
# na.action = function defining action to take for cases that have NA for
#     any of ftime, fstatus, cov1, cov2 cengroup, or subset.
# gtol = iteration stops when a function of the gradient is < gtol.
# maxiter = maximum # of iterations in Newton algorithm (0 computes
#     scores and var at init, but performs no iterations)
# init = initial values of regression parameters
###### data manipulation
  d = data.frame(ftime=ftime,fstatus=fstatus,
         strata=if (missing(strata)) 1 else strata)
# add the covariates to the data frame
 if (!missing(cov1)) {
    d$cov1 = as.matrix(cov1)
    nc1 = ncol(d$cov1)
  } else {nc1 = 0}
  if (!missing(cov2)) {
    d$cov2 = as.matrix(cov2)
    nc2 = ncol(d$cov2)
  } else {nc2 = 0}
    np = nc1 + nc2
# subset the data if specified
  if (!missing(subsets)) d = d[subsets,]
# obtain the number of cases with missing values
  tmp = nrow(d)
  d = na.action(d)
  if (nrow(d)!= tmp) cat(format(tmp-nrow(d)),
     'cases omitted due to missing values\n')
# generate some variables
  d = d[order(d$ftime),]
  usl = unique(d$strata)
  ns = length(usl)
  d$cenind = ifelse(d$fstatus==cencode,1,0)
  d$fstatus = ifelse(d$fstatus==failcode,1,2*(1-d$cenind))
  uft = sort(unique(d$ftime[d$fstatus==1]))
  ndf = length(uft)

  if (nc2 == 0) d$tfs = 0
  else d$tfs = tf(d$ftime)
 
# want censoring dist km at ftime
    if (ctype == 1) {
       d$subject = 1:length(d$ftime)
       gout = NULL
       for (j in 1:ns)  {
         dsub = subset(d, strata == usl[j])
         u =  do.call('survfit',list(formula=Surv(dsub$ftime,dsub$cenind) ~ 1))
         u = summary(u,times=sort(dsub$ftime*(1-.Machine$double.eps)))
         tmp = cbind(uuu=u$surv, subject=dsub$subject)
         gout = rbind(gout, tmp)  
       }
       d = merge(d,gout,by="subject")
    } else {
      u = do.call('survfit',list(formula= Surv(d$ftime,d$cenind) ~ 1))
      u = summary(u,times= sort(d$ftime* (1- .Machine$double.eps)))
      d$uuu = u$surv
   }

######## start of newton-raphson
  if (missing(init)) b = rep(0,np)
  else b = init
  stepf = .5

  for (ll in 0:maxiter) {
  
      ## function to calculate loglik, score, and info
      z  = crrfsvs(d=d,b=b,nc2=nc2,np=np,ns=ns,usl=usl)
  
      if (max(abs(z[[2]])* pmax(abs(b),1)) < max(abs(z[[1]]),1)* gtol) {
        converge = TRUE
        break
      }
  
      if (ll==maxiter) {
        converge = FALSE
        break
      }
  
      h = z[[3]]
      dim(h) = c(np,np)
      sc = solve(h,z[[2]])
      bn = b + sc
  
      ## function to calculate loglik
      fbn = crrfs(d=d,bn=bn,nc2=nc2,np=np,ns=ns,usl=usl)
  
      ## backtracking loop
      i = 0
      while (is.na(fbn) || fbn< z[[1]]+(1e-4)*sum(sc*z[[2]])) {
        i = i+ 1
        sc = sc* stepf
        bn = b+ sc
  
      ## function to calculate loglik
      fbn = crrfs(d=d,bn=bn,nc2=nc2,np=np,ns=ns,usl=usl)
        if (i>20) break
      }
  
      if (i>20) {
        converge = FALSE
        break
      }
  
      b = c(bn)
  }
# end of newton-raphson

### function to estimate variance (two matrices)
  if (ctype==1) v = crrvvs(d=d,b=b,nc2=nc2,np=np,ns=ns,usl=usl)
  if (ctype==2) v = .C("crrvvh", as.double(d$ftime),as.integer(d$fstatus),
          as.integer(length(d$ftime)), as.double(d$cov1), as.integer(np-nc2),
          as.double(d$cov2), as.integer(nc2), as.double(d$tfs), as.integer(ndf),
          as.integer(d$strata),as.integer(ns), as.double(d$uuu), as.double(b), 
          double(np*np), double(np*np),PACKAGE="crrSC")[14:15]
  dim(v[[2]]) = dim(v[[1]]) = c(np,np)
  h = solve(v[[1]])  
  v = h %*% v[[2]] %*% t(h)

# calculate the unique failure time and time varying covariates across strata
  uft = unique(ftime)
  tfs = 0
  if (nc2 != 0) tfs = tf(uft)

  z = list(coef=b, loglik=-z[[1]], score=-z[[2]], inf=matrix(z[[3]],np,np),
        var=v, h=h, uftime=uft, tfs=tfs, converged=converge)
  class(z) = 'crrs'
  z
}

################################################################################
# similar to "crrfsv" in "cmprsk" package (added stratification)
# obtain loglik, score, and information
crrfsvs = function(d=d,b=b,nc2=nc2,np=np,ns=ns,usl=usl)
{
  z1 = z2 = z3 = 0
  for (j in 1:ns)
  {
    # subsetting data 
    dsub = subset(d,strata==usl[j])
    ftime = dsub$ftime
    cenind = dsub$cenind
    fstatus = dsub$fstatus
    uft = sort(unique(ftime[fstatus==1]))
    ndf = length(uft)
    cov1 = dsub$cov1
    cov2 = dsub$cov2
    tfs = unique(subset(dsub$tfs,fstatus==1))
    uuu = dsub$uuu

    # calculation
    temp = .C("crrfsv", as.double(ftime),as.integer(fstatus),
      as.integer(length(ftime)), as.double(cov1), as.integer(np-nc2),
      as.double(cov2), as.integer(nc2), as.double(tfs), as.integer(ndf),
      as.double(uuu), as.double(b), double(1), double(np), double(np*np),PACKAGE="crrSC")[12:14]
    z1 = z1 + temp[[1]]
    z2 = z2 + temp[[2]]
    z3 = z3 + temp[[3]]
  }
  list(z1,z2,z3)
}


################################################################################

# similar to "crrf" in "cmprsk" package (added stratification)
# to obtain the loglik

crrfs = function(d=d,bn=bn,nc2=nc2,np=np,ns=ns,usl=usl)
{
  z4 = 0
  for (j in 1:ns)
  {
    # subset data and calculate
    dsub = subset(d,strata==usl[j])
    ftime = dsub$ftime
    cenind = dsub$cenind
    fstatus = dsub$fstatus
    uft = sort(unique(ftime[fstatus==1]))
    ndf = length(uft)
    cov1 = dsub$cov1
    cov2 = dsub$cov2
    tfs = unique(subset(dsub$tfs,fstatus==1))
    uuu = dsub$uuu

    # calculation
    temp = .C("crrf", as.double(ftime),as.integer(fstatus),
        as.integer(length(ftime)), as.double(cov1), as.integer(np-nc2),
        as.double(cov2), as.integer(nc2), as.double(tfs), as.integer(ndf),
        as.double(uuu), as.double(bn), double(1),PACKAGE="crrSC")[[12]]
     z4 = z4 + temp
  }
  z4
}


################################################################################
# similar to "crrvv" in "cmprsk" package (added stratification)
# to obtain the robust var estimate for regular stratified data

crrvvs = function(d = d,b = b,nc2=nc2,np=np,ns=ns,usl=usl)
{
  v1 = v2 = 0
  for (j in 1:ns)
  {
    # subset data and calculate
    dsub = subset(d,strata==usl[j])
    ftime = dsub$ftime
    cenind = dsub$cenind
    fstatus = dsub$fstatus
    uft = sort(unique(ftime[fstatus==1]))
    ndf = length(uft)
    cov1 = dsub$cov1
    cov2 = dsub$cov2
    tfs = unique(subset(dsub$tfs,fstatus==1))
    uuu = dsub$uuu

    # calculation
    temp = .C("crrvv", as.double(ftime),as.integer(fstatus),
        as.integer(length(ftime)), as.double(cov1), as.integer(np-nc2),
        as.double(cov2), as.integer(nc2), as.double(tfs), as.integer(ndf),
        as.double(uuu), as.double(b), double(np*np), double(np*np),PACKAGE="crrSC")[12:13]
    v1 = v1 + temp[[1]]
    v2 = v2 + temp[[2]]
  }
  list(v1,v2)
}                                                             

################################################################################
print.crrs = function (x, ...) 
{
    cat("convergence: ", x$converged, "\n")
    cat("coefficients:\n")
    print(signif(x$coef, 4), ...)
    v <- sqrt(diag(x$var))
    cat("standard errors:\n")
    print(signif(v, 4), ...)
    v <- 2 * (1 - pnorm(abs(x$coef)/v))
    cat("two-sided p-values:\n")
    print(signif(v, 2), ...)
    invisible()
}



