SparamW.S <-
function(X,y,delta,N,q,sigma0,MAXIT,TOL,ialg=3,seed=153) {
X <- as.matrix(X); n <- length(y); p <- ncol(X); b <- 0.5
#
# generate matrix of Beta-s
#
zbet  <- BtamatW(X,y,delta,N,q,MAXIT,TOL,seed=seed) 
beta  <- zbet$beta
iok      <- (1:N)[!is.na(beta[,1])]
beta     <- beta[iok,,drop=FALSE]
N        <- length(iok)
gamok <- 0
s0    <- 1; c0 <- -0.1351788
kappa <- 9.E9
smink <- 9.E9
mes2  <- integer(4)
  cat(".")
  flush.console()
  now <- Sys.time()
  now <- now+6
for (j in 1:N) {
 betj  <- matrix(beta[j,],nrow=N,ncol=p,byrow=TRUE)
 if (gamok==0) {
  gamj <- sweep(beta,2,beta[j,])
  nu   <- apply(gamj,1,Nrm2)}
 betgamj <- betj
 B      <- (1:N)[nu <= kappa]
 ni     <- length(B)
 if (ni==0) next
 betgamb <- betgamj[B,]
#
# initial scale s_n(betaj+gama)
#
 tmp    <- s.eq.logW (X,y,ni,delta,sigma0,b,betgamb,MAXIT,TOL,mu0=c0,s0=s0,
           ipsi=4,xk=1.717816,ialg=ialg,meth=4)
 sj     <- tmp$S
 mes2   <- mes2+tmp$mes2
#
# compute S_n(beta,gamma) and minimize over gamma
#
 Sj     <- rep(100000,N) 
 gamjb  <- gamj[B,]
 betjb  <- betj[B,]
 tmp    <- S.eq.logW (X,y,ni,delta,sj,sigma0,b,betjb,gamjb,MAXIT,TOL,mu0=c0,
           ipsi=4,xk=1.717816,ialg=ialg)
 Sj[B]  <- tmp$S
 mes2   <- mes2+tmp$mes2
 omega  <- min(Sj)
 kStar  <- (1:N)[Sj==omega] 
 numin  <- nu[kStar]
 nustar <- min(numin)
 kstar  <- (kStar)[numin==nustar]
# cat(j,kstar,round(c(beta[j,2],gama[kstar,2],omega,nustar,kappa),4))
#
# Check for equal norm
#
 if (abs(nustar-kappa)<1e-6 & omega>smink) next
 Bbar   <- setdiff(1:N,B)
 truemin <- 1
 if (length(Bbar)>0) for (k in Bbar) {
  tmp <- s.eq.logW (X,y,1,delta,sigma0,b,betgamj[k,],MAXIT,TOL,mu0=c0,s0=s0,
         ipsi=4,xk=1.717816,ialg=ialg,meth=4)
  sk  <- tmp$S; mes2 <- mes2+tmp$mes2
  tmp <- S.eq.logW (X,y,1,delta,sk,sigma0,b,betj[k,],gamj[k,],MAXIT,TOL,mu0=c0,
         ipsi=4,xk=1.717816,ialg=ialg)
  tst <- tmp$S; mes2 <- mes2+tmp$mes2
  if (tst < omega) {truemin <- 0; break}
 }
 if (truemin==1) {
  smink  <- omega
  Jstar  <- j
  Kstar  <- kstar
  ind    <- match(kstar,B,nomatch=0)
  Smin   <- sj[ind]
  smin   <- Sj[kstar]
  gamin  <- gamj
  kappa  <- nustar} 
    if (Sys.time() >= now) {cat("."); flush.console(); now <- now+6}
} # end for j
#
  cat("\n"); flush.console()
cmes2  <- c(as.character(mes2)," ")
Smin   <- Smin[1]
Eqexit <- paste(c("Normal (","Nit=Maxit (","f(a)*f(b)>0 (","|f(a)-f(b)|<tl ("," "),
          cmes2,sep="",collapse=")  ")
Bmin0  <- beta[Jstar,]
Smin0  <- Smin
zres   <- list(Bmin=Bmin0,Smin=Smin0,smin=smin,beta=beta,gama=gamin,kappa=kappa,
          Jstar=Jstar,Kstar=Kstar,Eqexit=Eqexit)
zres}

