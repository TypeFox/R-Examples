ds <-
function(z, Fm=6, Fs=6, type=c("daily", "monthly"), searchQ=TRUE, lag.max=20, ic=c("BIC","AIC"),
	standardizeQ=TRUE){
    if (is.matrix(z)&&ncol(z)>1) stop("Error z: must be univariate time series, vector or matrix with one column")
    if (is.matrix(z)) z <- as.vector(z)
    if (is.ts(z)&&frequency(z)==12) type <- "monthly" else type <- match.arg(type)
    if (type=="daily") s<-365.25 else s<-12
    if (! (s%in%c(12,365.25))) stop("Error: s must be 12 or 365.25")
    if (min(Fm,Fs) < 0) stop("Error: both Fm and Fs must be >= 0")
    if (max(Fm,Fs) > 6) stop("Error: Fm and Fs are restricted to be <= 6")
    if (max(Fm,Fs) > s/2) stop("error: Fm or Fs setting too large")

    ic <- match.arg(ic)
    nz <- length(z)
    if (nz < s) stop("error: need at least one year (366 days or 12 months)")
#
if (!searchQ) { #only one specific model
  ans <- getds(z=z, s=s, Fm=Fm, Fs=Fs, ic=ic, lag.max=lag.max, standardizeQ=standardizeQ)
  m <- matrix(ans$dspar, nrow=1, ncol=4)
  zds <- ans$z
  } else {
#enumerate all models, find best
      maxFm <- Fm
      maxFs <- Fs
      m <- matrix(numeric(4*((maxFm+1)*(maxFs+1))), ncol=4)
      colnames(m) <- c("Fm", "Fs", "p", ic)
      i <- 0
      for (iFm in 0:maxFm)
        for (iFs in 0:maxFs) {
          i <- i+1
          m[i,] <- getds(z=z, s=s, Fm=iFm, Fs=iFs, ic=ic, lag.max=lag.max, standardizeQ=standardizeQ)$dspar
          }
      }
      rownames(m) <- rep(" ", nrow(m))
      ind<-which.min(m[,4])
      rownames(m)[ind]<-"*" #best one
      FmOpt <- m[ind,1]
      FsOpt <- m[ind,2]
      zds <- getds(z=z, s=s, Fm=FmOpt, Fs=FsOpt, ic=ic, lag.max=lag.max, standardizeQ=standardizeQ)$z
mnames <- colnames(m)
out <- list(dspar=m, z=zds)
class(out) <- "deseasonalize"
out
}
