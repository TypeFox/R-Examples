rivGibbs=function(Data,Prior,Mcmc) {
#
# revision history:
#    R. McCulloch original version 2/05 
#    p. rossi 3/05 
#    p. rossi 1/06 -- fixed error in nins
#    p. rossi 1/06 -- fixed def Prior settings for nu,V
#    3/07 added classes
#    W. Taylor 4/15 - added nprint option to MCMC argument
#
#
# purpose: 
#   draw from posterior for linear I.V. model
#
# Arguments:
#   Data -- list of z,w,x,y
#        y is vector of obs on lhs var in structural equation
#        x is "endogenous" var in structural eqn
#        w is matrix of obs on "exogenous" vars in the structural eqn
#        z is matrix of obs on instruments
#   Prior -- list of md,Ad,mbg,Abg,nu,V
#        md is prior mean of delta
#        Ad is prior prec
#        mbg is prior mean vector for beta,gamma
#        Abg is prior prec of same
#        nu,V parms for IW on Sigma
#
#   Mcmc -- list of R,keep 
#        R is number of draws
#        keep is thinning parameter
#        nprint - print estimated time remaining on every nprint'th draw
#
#   Output: 
#      list of draws of delta,beta,gamma and Sigma
# 
#   Model:
#
#    x=z'delta + e1
#    y=beta*x + w'gamma + e2
#        e1,e2 ~ N(0,Sigma)
#
#   Priors
#   delta ~ N(md,Ad^-1)
#   vec(beta,gamma) ~ N(mbg,Abg^-1)
#   Sigma ~ IW(nu,V)
#
#   check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of z,w,x,y")}
    if(is.null(Data$z)) {pandterm("Requires Data element z")}
    z=Data$z
    if(is.null(Data$w)) {pandterm("Requires Data element w")}
    w=Data$w
    if(is.null(Data$x)) {pandterm("Requires Data element x")}
    x=Data$x
    if(is.null(Data$y)) {pandterm("Requires Data element y")}
    y=Data$y
#
# check data for validity
#
if(!is.vector(x)) {pandterm("x must be a vector")}
if(!is.vector(y)) {pandterm("y must be a vector")}
n=length(y)
if(!is.matrix(w)) {pandterm("w is not a matrix")}
if(!is.matrix(z)) {pandterm("z is not a matrix")}
dimd=ncol(z)
dimg=ncol(w)
if(n != length(x) ) {pandterm("length(y) ne length(x)")}
if(n != nrow(w) ) {pandterm("length(y) ne nrow(w)")}
if(n != nrow(z) ) {pandterm("length(y) ne nrow(z)")}
#
# check for Prior
#
if(missing(Prior))
   { md=c(rep(0,dimd));Ad=BayesmConstant.A*diag(dimd); 
     mbg=c(rep(0,(1+dimg))); Abg=BayesmConstant.A*diag((1+dimg));
     nu=3; V=diag(2)}
else
   {
    if(is.null(Prior$md)) {md=c(rep(0,dimd))} 
       else {md=Prior$md}
    if(is.null(Prior$Ad)) {Ad=BayesmConstant.A*diag(dimd)} 
       else {Ad=Prior$Ad}
    if(is.null(Prior$mbg)) {mbg=c(rep(0,(1+dimg)))} 
       else {mbg=Prior$mbg}
    if(is.null(Prior$Abg)) {Abg=BayesmConstant.A*diag((1+dimg))} 
       else {Abg=Prior$Abg}
    if(is.null(Prior$nu)) {nu=3}
       else {nu=Prior$nu}
    if(is.null(Prior$V)) {V=nu*diag(2)}
       else {V=Prior$V}
   }
#
# check dimensions of Priors
#
if(ncol(Ad) != nrow(Ad) || ncol(Ad) != dimd || nrow(Ad) != dimd) 
   {pandterm(paste("bad dimensions for Ad",dim(Ad)))}
if(length(md) != dimd)
   {pandterm(paste("md wrong length, length= ",length(md)))}
if(ncol(Abg) != nrow(Abg) || ncol(Abg) != (1+dimg) || nrow(Abg) != (1+dimg)) 
   {pandterm(paste("bad dimensions for Abg",dim(Abg)))}
if(length(mbg) != (1+dimg))
   {pandterm(paste("mbg wrong length, length= ",length(mbg)))}
#
# check MCMC argument
#
if(missing(Mcmc)) {pandterm("requires Mcmc argument")}
else
   {
    if(is.null(Mcmc$R)) 
       {pandterm("requires Mcmc element R")} else {R=Mcmc$R}
    if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
    if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
      if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
   }

#
# print out model
#
cat(" ",fill=TRUE)
cat("Starting Gibbs Sampler for Linear IV Model",fill=TRUE)
cat(" ",fill=TRUE)
cat(" nobs= ",n,"; ",ncol(z)," instruments; ",ncol(w)," included exog vars",fill=TRUE)
cat("     Note: the numbers above include intercepts if in z or w",fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("mean of delta ",fill=TRUE)
print(md)
cat("Adelta",fill=TRUE)
print(Ad)
cat("mean of beta/gamma",fill=TRUE)
print(mbg)
cat("Abeta/gamma",fill=TRUE)
print(Abg)
cat("Sigma Prior Parms",fill=TRUE)
cat("nu= ",nu," V=",fill=TRUE)
print(V)
cat(" ",fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat(" ",fill=TRUE)

###################################################################
# Keunwoo Kim
# 09/03/2014
###################################################################
draws=rivGibbs_rcpp_loop(y, x, z, w, mbg, Abg, md, Ad, V, nu, R, keep, nprint)
###################################################################

attributes(draws$deltadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$deltadraw)$mcpar=c(1,R,keep)
attributes(draws$betadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$betadraw)$mcpar=c(1,R,keep)
attributes(draws$gammadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$gammadraw)$mcpar=c(1,R,keep)
attributes(draws$Sigmadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(draws$Sigmadraw)$mcpar=c(1,R,keep)

return(draws)
}
