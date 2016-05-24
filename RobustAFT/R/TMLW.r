TMLW <-
function(X,y,delta,Beta.t,sigma.t,const,tl,tu,control) {
# TML, logWeibull case, iteratively reweighting algorithm
maxit.sigma   <- control$maxit.sigma
tol.sigma     <- control$tol.sigma
maxit.Beta    <- control$maxit.Beta
tol.Beta      <- control$tol.Beta
Maxit.TML     <- control$Maxit.TML
tol.TML.sigma <- control$tol.TML.sigma
tol.TML.Beta  <- control$tol.TML.Beta
alg.sigma     <- control$alg.sigma
nitmon        <- control$nitmon
p <- ncol(X)
Nit    <- 1
Mxt <- 10
Beta0  <- Beta.t
sigma0 <- sigma.t
mui.t  <- X %*% as.matrix(Beta.t)
rs0    <- (y-mui.t)/sigma.t
wgt    <- ww(rs0,tl, tu)
DeltaB <- Deltas <- 100000
Nit    <- 1
while ( ( Deltas > tol.TML.sigma | DeltaB > tol.TML.Beta ) & Nit <= Maxit.TML) {
mui   <- X %*% as.matrix(Beta0)                                                            # check
if (alg.sigma ==1) {
  zs    <- TML.SigmaW(X,y,delta,sigma0,sigma.t,mui,mui.t,wgt,const,tl,tu,tol=tol.sigma,maxit=maxit.sigma,nitmon)
  sigma <- zs$sigma; nits <- zs$nit }
if (alg.sigma !=1) {
# determine lower and upper values for sigma
  smu  <- sml  <- sigma0;  nit  <- 1 
  funu <- funl <- TML.Ave2W(X,y,delta, smu,sigma.t, mui,mui.t, wgt,tl,tu)-const
  while(funu > 0 & nit <= Mxt) {
    smu  <- smu*1.5
    funu <-       TML.Ave2W(X,y,delta, smu,sigma.t, mui,mui.t, wgt,tl,tu)-const
    nit  <- nit+1}
  if(nitmon) cat("smu, funu, nit",smu,funu,nit,"\n")
  nit <- 1
  while(funl < 0 & nit <= Mxt) {
    sml  <- sml*0.5
    funl <-       TML.Ave2W(X,y,delta, sml,sigma.t, mui,mui.t, wgt,tl,tu)-const
    nit  <- nit+1}
  if(nitmon) cat("sml, funl, nit",sml,funl,nit,"\n")
# compute new sigma
  zs <- regfal(TML.Ave2W,cc=const,lower=sml,upper=smu,nint=5,tol=tol.sigma,maxit=maxit.sigma,
               X=X,y=y,delta=delta, sigma.t=sigma.t, mui=mui, mui.t=mui.t, wgt=wgt,cl=tl,cu=tu)
  sigma  <- zs$solution; nits <- zs$nit}
  if(nitmon) cat("sigma", nits, "\n")
# compute new Beta
  zB   <- TML.BetaW(X,y,delta,Beta0,sigma,Beta.t,sigma.t, tl,tu,maxit=maxit.Beta,tol=tol.Beta,nitmon)
  Beta <- zB$Beta
  DeltaB <- max(abs(Beta - Beta0)); Deltas <- abs(sigma-sigma0)
  if(nitmon) cat("Nit=",Nit,"DeltaB",DeltaB,"Deltas",Deltas, "\n")
  Beta0 <- Beta; sigma0 <- sigma
  Nit   <- Nit+1}
zres <- list(Beta=Beta,sigma=sigma,Nit=Nit)
zres}

