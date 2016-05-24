RefSG <-
function(X,y,delta,Bmin0,Smin0,control) {
# Refinement for parametric S, Gaussian case, iteratively reweighting algorithm
maxit.sigma <- control$maxit.sigma
tol.sigma   <- control$tol.sigma
maxit.Beta  <- control$maxit.Beta
tol.Beta    <- control$tol.Beta
Maxit.SG    <- control$Maxit.S
tol.SG.sigma<- control$tol.S.sigma
tol.SG.Beta <- control$tol.S.Beta
alg.sigma   <- control$alg.sigma
nitmon      <- control$nitmon
p <- ncol(X);  Nit    <- 1; Mxt <- 10
Bmin   <- Bmin0; Smin <- Smin0
DeltaB <- DeltaS <- 100000
while ( ( DeltaS > tol.SG.sigma | DeltaB > tol.SG.Beta ) & Nit <= Maxit.SG) {
if (alg.sigma ==1) {
  zS <- RefSigmaG(Smin,Bmin,X,y,delta,tol=tol.sigma,maxit=maxit.sigma,nitmon)
  Smin <- zS$sigma; nitS <- zS$nit }
if (alg.sigma !=1) {
# determine lower and upper values for sigma
  smu  <- sml  <- Smin;  nit  <- 1 
  funu <- funl <- RefAve2G(smu,Bmin,X,y,delta)-0.5
  while(funu > 0 & nit <= Mxt) {
    smu  <- smu*1.5
    funu <- RefAve2G(smu,Bmin,X,y,delta)-0.5
    nit  <- nit+1}
  if (nitmon) cat("smu, funu, nit",smu,funu,nit,"\n")
  nit <- 1
  while(funl < 0 & nit <= Mxt) {
    sml  <- sml*0.5
    funl <- RefAve2G(sml,Bmin,X,y,delta)-0.5
    nit  <- nit+1
  }
  if (nitmon) cat("sml, funl, nit",sml,funl,nit,"\n")
# compute new Smin
  zS <- regfal(RefAve2G,cc=0.5,lower=sml,upper=smu,nint=5,tol=tol.sigma,maxit=maxit.sigma,
               Beta=Bmin,X=X,y=y,delta=delta)
  Smin  <- zS$solution; nitS <- zS$nit}
  if (nitmon) cat("Smin", Smin, nitS, "\n")
# compute new Bmin
  zB   <- RefBetaG(X,y,delta,Bmin,Smin,maxit=maxit.Beta,tol=tol.Beta,nitmon)
  Bmin <- zB$Beta
  DeltaB <- max(abs(Bmin - Bmin0)); DeltaS <- abs(Smin-Smin0)
  if (nitmon) cat("Nit=",Nit,"DeltaB",DeltaB,"DeltaS",DeltaS,"\n")
  Bmin0 <- Bmin; Smin0 <- Smin
  Nit   <- Nit+1}
zres <- list(Bmin=Bmin,Smin=Smin,Nit=Nit)
zres}

