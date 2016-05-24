# Functions to test properties of random proportions

# ************************** Utilities *******************************
# Composition operator
# Converts vector v to proportions

local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

comp <- function(v) {
  v <- v[v>0 & !is.na(v)]
  y <- abs(v)/sum(abs(v)); y; };

panel.hist <- function(x, ...) {
   usr <- par("usr"); on.exit(par(usr))
   h <- hist(x, breaks="Sturges", plot=FALSE)
   breaks <- h$breaks; nB <- length(breaks)
   y <- h$counts; y <- y/sum(y)
   par(usr = c(usr[1:2], 0, max(y)*1.5) )
   rect(breaks[-nB], 0, breaks[-1], y, col=switch(getWinVal()$MDL,"lightblue1","mistyrose","darkseagreen1"))
   box() }

# ************************** Multinomial *****************************

# Multinomial random sample
#   Input:  ns = number of sample
#           N = sample size
#           pvec = probability vector (sums to 1)
#           prop = report format:
#                  counts (prop=F) or proportions (prop=T)
#   Output: matrix (ns x length(pvec))
#           each row is a multinomial sample that sums to N

rmul <- function(ns,N,pvec,prop=T) {
  pvec <- abs(pvec)/sum(abs(pvec)); # forces sum to 1
  np   <- length(pvec);
  ymat <- rmultinom(ns,N,pvec);
  yout <- t(ymat);
  if(prop) yout <- sweep( yout, 1, N, "/");
  yout; };

# ***** Test of rmul *****

testrmul <- function(prop=T,ndec=4) {
  getWinVal(scope="L");
  pvec <- comp(pvec); np <- length(pvec);
  y <- rmul(ns,N,pvec,prop);
  dimnames(y)[[2]] <- paste("p",1:np,sep="");
  ym <- apply(y,2,mean);
  yv <- apply(y,2,var);
  yvv <- pvec*(1-pvec)/N;
  ys <- apply(y,2,sd);
  PVEC <- show0(round(c(pvec,rep(0,6-np)),ndec),ndec);
  YM   <- show0(round(c(ym,rep(0,6-np)),ndec),ndec);
  YS <- show0(round(c(ys,rep(0,6-np)),ndec),ndec);
  setWinVal(list(pvec=PVEC,ym=YM,ys=YS));
  pairs(y,pch=16,cex=0.5,col="blue",gap=0,diag.panel=panel.hist);
  writeList(list(ns=ns,N=N,pvec=pvec,y=y),fnam="RpData.txt",format="P");
  invisible(y); };

# ************************** Dirichlet *******************************

# Dirichlet random sample
#   Input:  ns = number of samples
#           N = effective sample size
#           pvec = probability vector (sums to 1)
#   Output: matrix (ns x length(pvec))
#           each row is a Dirichlet sample that sums to 1

rdir <- function(ns,N,pvec) {
  pvec <- comp(pvec); # forces sum to 1
  np <- length(pvec);
  yvec <- rgamma(ns*np,shape=N*pvec,scale=N);
  y1 <- matrix(yvec,nrow=ns,ncol=np,byrow=T);
  ys <- apply(y1,1,"sum");
  y2 <- sweep(y1,1,ys,"/");
  y2; };

# ***** Test of rdir *****

testrdir <- function(ndec=4) {
  getWinVal(scope="L");
  pvec <- comp(pvec); np <- length(pvec);
  y   <- rdir(ns,N,pvec); 
  dimnames(y)[[2]] <- paste("p",1:np,sep="");
  ym  <- apply(y,2,mean);
  yv  <- apply(y,2,var);
  yvv <- pvec*(1-pvec)/(N+1);
  ys  <- apply(y,2,sd);
  PVEC <- show0(round(c(pvec,rep(0,6-np)),ndec),ndec);
  YM   <- show0(round(c(ym,rep(0,6-np)),ndec),ndec); 
  YS <- show0(round(c(ys,rep(0,6-np)),ndec),ndec);
  setWinVal(list(pvec=PVEC,ym=YM,ys=YS));
  pairs(y,pch=16,cex=0.5,col="red",gap=0,diag.panel=panel.hist);
  writeList(list(ns=ns,N=N,pvec=pvec,y=y),fnam="RpData.txt",format="P");
  invisible(y); };

# ************************ Logit-normal ******************************

# Logit-normal random sample
#   Input:  ns = number of samples
#           sig = standard error
#           pvec = probability vector (sums to 1)
#   Output: matrix (ns x length(pvec))
#           each row is a logit-normal sample that sums to 1

rlgtnorm <- function(ns,sig,pvec) {
  pvec <- comp(pvec);    # forces sum to 1
  #pvec <- pvec[pvec>0];  # omits zero elements
  np   <- length(pvec);
  mvec <- log(pvec);
  vvec <- rnorm(ns*np,mean=mvec,sd=sig);
  y1   <- matrix(exp(vvec),nrow=ns,ncol=np,byrow=T);
  ys   <- apply(y1,1,"sum");
  y2   <- sweep(y1,1,ys,"/");
  y2; };

# Test of rlgtnorm

testrln <- function(ndec=4) {
  getWinVal(scope="L");
  pvec <- comp(pvec); np <- length(pvec);
  y   <- rlgtnorm(ns,sig,pvec);
  dimnames(y)[[2]] <- paste("p",1:np,sep="");
  ym  <- apply(y,2,mean);
  yv  <- apply(y,2,var);
  yvv <- sig^2*pvec^2*(1 - 2*pvec + sum(pvec^2));
  ys  <- apply(y,2,sd);
  iy  <- sort( sample(1:np, 2, replace=F) );
  i1  <- iy[1]; i2 <- iy[2];
  y1  <- y[,i1]; y2 <- y[,i2];
  y12m <- mean(log(y1/y2)); y12 <- log(pvec[i1]/pvec[i2]);
  v12m <- var(log(y1/y2)); v12 <- 2*sig^2;
  PVEC <- show0(round(c(pvec,rep(0,6-np)),ndec),ndec);
  YM   <- show0(round(c(ym,rep(0,6-np)),ndec),ndec);
  YS <- show0(round(c(ys,rep(0,6-np)),ndec),ndec);
  setWinVal(list(pvec=PVEC,ym=YM,ys=YS));
  pairs(y,pch=16,cex=0.5,col="forestgreen",gap=0,diag.panel=panel.hist);
  writeList(list(ns=ns,sig=sig,pvec=pvec,y=y),fnam="RpData.txt",format="P");
  invisible(y); };

getMod <- function(){
  getWinVal(scope="L");
  rmod <- switch(MDL,"testrmul","testrdir","testrln")
  get(rmod)()
}
if (!require(PBSmodelling, quietly=TRUE)) stop("The PBSmodelling package is required for this example")
createWin("RanPropWin.txt")

}) # end local scope

