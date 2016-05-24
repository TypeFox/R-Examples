local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

simData <- function() {

  # calculate parameters
  getWinVal(scope="L"); Tcv <- Tsd/Tmn; 
  esig2 <- 1 + Tcv^2; Lsig <- sqrt(log(esig2));
  Lmu <- log(Tmn/sqrt(esig2));
  Gshp <- 1/Tcv^2; Gscl <- Tmn/Gshp;
  Npar <- list(Tcv=Tcv,Nmn=Tmn,Nsd=Tsd);
  setWinVal(lapply(Npar,round,4));
  Lpar <- list(Lmn=Tmn,Lsd=Tsd,Lmu=Lmu,Lsig=Lsig);
  setWinVal(lapply(Lpar,round,4));
  Gpar <- list(Gmn=Tmn,Gsd=Tsd,Gshp=Gshp,Gscl=Gscl);
  setWinVal(lapply(Gpar,round,4));

  # generate samples
  RVx <- rnorm(ns,Tmn,Tsd); tput(RVx)
  RVy <- rlnorm(ns,Lmu,Lsig); tput(RVy)
  RVz <- rgamma(ns,shape=Gshp,scale=Gscl); tput(RVz)
  
  # calculate estimates
  SNpar <- list(SNmn=mean(RVx),SNsd=sd(RVx));
  setWinVal(lapply(SNpar,round,4));
  SLmn <- mean(RVy); SLsd <- sd(RVy);
  logy <- log(RVy); SLmu <- mean(logy); SLsig <- sd(logy);
  SLpar <- list(SLmn=SLmn,SLsd=SLsd,SLmu=SLmu,SLsig=SLsig);
  setWinVal(lapply(SLpar,round,4));
  SGmn <- mean(RVz); SGsd <- sd(RVz); SGcv <- SGsd/SGmn;
  SGshp <- 1/SGcv^2; SGscl <- SGmn/Gshp;
  SGpar <- list(SGmn=SGmn,SGsd=SGsd,SGshp=SGshp,SGscl=SGscl);
  setWinVal(lapply(SGpar,round,4));
  invisible(); }    
  
plDens <- function() {
  resetGraph(); getWinVal(scope="L");
  qq <- c(0.001,0.999);
  xrng <- range( qnorm(qq,Tmn,Tsd), qlnorm(qq,Lmu,Lsig),
                 qgamma(qq,shape=Gshp,scale=Gscl)); 
  x1 <- seq(xrng[1],xrng[2],length=500);
  x2 <- seq(1e-04,xrng[2],length=500);
  y1 <- dnorm(x1,Tmn,Tsd);
  y2 <- dlnorm(x2,Lmu,Lsig);
  y3 <- dgamma(x2,shape=Gshp,scale=Gscl);
  yrng <- range(y1,y2,y3);
  yrng[2] <- min(yrng[2],2*max(y1,y2));
  plot(x1,y1,xlim=xrng,ylim=yrng,xlab="x",ylab="pdf",
    col="red",type="l",lwd=2);
  lines(x2,y2,col="forestgreen",lwd=2);
  lines(x2,y3,col="blue",lwd=2);
  
  legend("topright", legend = c("Normal", "Lognormal", "Gamma"),
         text.width = strwidth("1,000,000"), lwd=2, lty = 1,
         xjust = 1, yjust = 1, col = c("red","forestgreen","blue"))
  
  invisible(); }

plCum <- function() {
  resetGraph(); getWinVal(scope="L");
  qq <- seq(0.0001,0.9999,length=300);
  x <- qnorm(qq,Tmn,Tsd); y <- qlnorm(qq,Lmu,Lsig);
  z <- qgamma(qq,shape=Gshp,scale=Gscl);
  xrng <- range(c(x,y,z));
  par(mfrow=c(2,1))
  plot(x,qq,xlim=xrng,ylim=c(0,1),xlab="x",ylab="Q",
    col="red",type="l",lwd=2);
  lines(y,qq,col="forestgreen",lwd=2);
  lines(z,qq,col="blue",lwd=2);
  
  
  
  legend("bottomright", legend = c("Normal", "Lognormal", "Gamma"),
         text.width = strwidth("1,000,000"), lwd=2, lty = 1, title="Theoretical",
         xjust = 1, yjust = 1, col = c("red","forestgreen","blue"))
  
  tget(RVx); x1 <- sort(RVx);
  tget(RVy); x2 <- sort(RVy);
  tget(RVz); x3 <- sort(RVz);
  qq <- (1:ns)/ns;
  xrng <- range(c(x1,x2,x3));
  plot(x1,qq,xlim=xrng,ylim=c(0,1),xlab="x obs",ylab="Q",
    col="red",type="l",lwd=2);
  lines(x2,qq,col="forestgreen",lwd=2);
  lines(x3,qq,col="blue",lwd=2);

  legend("bottomright", legend = c("Normal", "Lognormal", "Gamma"),
         text.width = strwidth("1,000,000"), lwd=2, lty = 1, title="Simulated",
         xjust = 1, yjust = 1, col = c("red","forestgreen","blue"))

  invisible(); }

panel.hist <- function(x, ...) {
   usr <- par("usr"); on.exit(par(usr))
   h <- hist(x, breaks="Sturges", plot=FALSE)
   breaks <- h$breaks; nB <- length(breaks)
   y <- h$counts; y <- y/sum(y)
   par(usr = c(usr[1:2], 0, max(y)*1.5) )
   rect(breaks[-nB], 0, breaks[-1], y, col="red")
   box() }

plPair <- function() {
	tget(RVx); tget(RVy); tget(RVz); 
  pairs(list(Normal=RVx,Lognormal=RVy,Gamma=RVz),
    pch=16,cex=0.5,col="forestgreen",gap=0,
    diag.panel=panel.hist); }

if (!require(PBSmodelling, quietly=TRUE)) stop("The PBSmodelling package is required for this example")
createWin("RanVarsWin.txt")

}) # end local scope
