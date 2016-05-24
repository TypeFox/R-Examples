getds <- function(z, s, Fm=6, Fs=6, ic=c("BIC","AIC"), lag.max=20, standardizeQ=TRUE){
  if (! (s%in%c(12,365.25))) stop("Error: s must be 12 or 365.25")
  if (min(Fm,Fs) < 0) stop("Error: both Fm and Fs must be >= 0")
  if (max(Fm,Fs) > 6) stop("Error: Fm and Fs are restricted to be <= 6")
#
  ic <- match.arg(ic)
  ICBest <- AdmissPenalty <- Inf
  pHat <- zds <- NA
  nz <- length(z)
  xt <- 1:nz
  if (Fm == 0) {
    zdm <- z-mean(z)
    } else if (Fm > 0) {
    X <- matrix(rep(1,(2*Fm+1)*nz), ncol=2*Fm+1)
    jj <- 2
    for (j in 1:Fm){
      x <- j*(2*pi/s)*xt
      X[,jj] <- sin(x)
      X[,jj+1] <- cos(x)
      jj <- jj+2
      }
    zdm <- resid(lm.fit(x=X, y=z))
    }
  estsd <- rep(sd(zdm), nz)
  J2 <- 0
  #
  if (!standardizeQ) fitsds <- rep(1, nz) else {
    if (Fs == 0) fitsds <- rep(var(z), nz) else if (Fs > 0) {
	    zsq <- zdm^2
    	X <- matrix(rep(1,(2*Fs+1)*nz), ncol=2*Fs+1)
    	jj <- 2
    	for (j in 1:Fs){
      		x <- j*(2*pi/s)*xt
      		X[,jj] <- sin(x)
      		X[,jj+1] <- cos(x)
      		jj <- jj+2
      		}
    	fitsds <- fitted(lm.fit(x=X, y=zsq))
  	}
  }
  if (all(fitsds>0)) {
    AdmissPenalty<-0  #otherwise set to Inf
    estsd <- sqrt(fitsds)
    J2 <- -sum(log(estsd))
    zds <- zdm/estsd
    ans <- SelectModel(zds, lag.max = lag.max, Criterion = ic, Best=2)
    pHat <- ans[1,1]
    ICBest <- ans[1,2]
    }
  if (ic == "BIC") parPenalty <- log(nz)*2*(Fm+Fs) else
       parPenalty <- 4*(Fm+Fs)
  BestIC <- ICBest -2*J2 + parPenalty + AdmissPenalty
  out<-c(Fm, Fs, pHat, BestIC)
  if (ic == "BIC") names(out) <- c("Fm", "Fs", "p", "BIC") else
        names(out) <- c("Fm", "Fs", "p", "AIC")
  list(dspar=out, z=zds)
}