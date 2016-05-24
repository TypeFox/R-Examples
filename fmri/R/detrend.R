fmri.detrend <- function(data,degree=1,accoef=0) {
  if (!class(data) == "fmridata") {
    warning("fmri.lm: data not of class <fmridata>. Try to proceed but strange things may happen")
  }
  cat("Start trend removal \n")
  ttt <- extract.data(data)
  dimttt <- dim(ttt)
  mask <- data$mask
  if (length(dimttt) != 4) {
    stop("Hmmmm, this does not seem to be a fMRI time series. I better stop executing! Sorry!\n")
  }
  n <- dimttt[4]
  z <- rep(1,n)
  if(degree>0) z <- cbind(rep(1,n),poly(1:n,degree))
  u <- svd(z,nv=0)$u
  dim(ttt) <- c(prod(dimttt[1:3]),dimttt[4])
  ttt[mask,] <- ttt[mask,] - ttt[mask,]%*%u%*%t(u)
  dim(ttt) <- dimttt
  cat("Finished trend removal \n")
  if(accoef>0){
  cat("Start prewhitening \n")
     rho0 <- 1/sqrt(1-accoef^2)
     rho1 <- accoef*rho0
     dim(ttt) <- c(prod(dimttt[1:3]),dimttt[4])
     ttt[mask,-1] <- rho0*ttt[mask,-1] - rho1*ttt[mask,-n]
     dim(ttt) <- dimttt
  cat("Finished prewhitening  \n")
  }
  data$ttt <- writeBin(as.numeric(ttt),raw(),4)
  invisible(data)
}
