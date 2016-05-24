#' Function to compute seasonal wave.
#'
#' @name compwaveconv
#' @title Seasonal Wave Computation
#' @param cmaxt is the time of the maximum chemical concentration, decimal
#' time in years.
#' @param jmod is the choice of model or pulse input function, an integer 
#' 1 through 14.
#' @param hlife is the model half-life in months, 1 to 4 months
#' @param mclass has not been implemented yet, but will provide
#' additional model options.
#' @return a numeric vector of size 361 with discrete values of the 
#' seasonal wave for decimal season \code{seq(0,1,1/360)}.
#' @note The seasonal wave is a dimensionless, periodic function of time 
#' with an annual cycle, similar to a mixture of sine and cosine functions 
#' often used to model seasonality in concentration data. However, the 
#' seasonal wave is better suited for modeling seasonal behavior of 
#' pesticide data than a mixture of sines and cosines. The pulse input 
#' function, represented by jmod, has either one or two distinct 
#' application seasons (when pesticides may be transported to the stream) 
#' of lengths from 1 to 6 months.  Therefore, 56 (14x4) choices for the 
#' wave function are available.
#' The numeric vector is a discrete approximation of the continuous wave 
#' function defined on the interval 0 to 1.
#' @keywords datagen
#' @author Aldo V. Vecchia
#' @export
#' @references Vecchia, A.V., Martin, J.D. and Gilliom, R.J., 2008, 
#' Modeling variability and trends in pesticide concentrations in streams: 
#' JAWRA Journal of the American Water Resources Association, v. 44, p. 
#' 1308--1324, 
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1752-1688.2008.00225.x/abstract}.
#' @examples
#' # evaluate seasonal wave for specified decimal seasons
#' # these example decimal dates represent days at points 0.25, 0.5, and 
#' # 0.75 percent of the way through the year and the end of the year
#' dseas <- c(0.25, 0.5, 0.75, 1)
#' swave <- compwaveconv(cmaxt=0.483, jmod=2, hlife=4, mclass=1)
#' swave[floor(360 * dseas)]
#' plot(seq(0,1,1/360),swave, typ="l")
compwaveconv <- function(cmaxt, jmod, hlife, mclass=1) {
  del <- 1 / 12
  txx <- seq(0, 1, 1 / 360)
  
  # phi is the decay rate corresponding with an approximate model 
  # half-life of 12 divided by phi months.
  if(hlife==1)  {
    phi <- 12
  } else if(hlife==2) { 
    phi <- 6
  } else if(hlife==3) { 
    phi <- 4
  } else if(hlife==4) {
    phi <- 3 
  } else { 
      stop("Half life value of ", hlife, 
           " is invalid, must be an integer 1 to 4.")
  }
  
  # wtx is the pulse input function and is greater than 0 during 
  # specified application sesons(s) and 0 otherwise.
  if(jmod==1) {
    wtx <- c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
    pkt <- 6 / 12
  } else if(jmod==2) {
    wtx <- c(0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0)
    pkt <- 7 / 12
  } else if(jmod==3) {
    wtx <- c(0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0)
    pkt <- 8 / 12
  } else if(jmod==4) {
    wtx <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
    pkt <- 9 / 12
  } else if(jmod==5) { 
    wtx <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0)
    pkt <- 9 / 12
  } else if(jmod==6) {
    wtx <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0)
    pkt <- 9 / 12
  } else if(jmod==7) {
    wtx <- c(0.5, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0.5)
    pkt <- 5 / 12
  } else if(jmod==8) {
    wtx <- c(0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0.5, 0.5)
    pkt <- 5 / 12
  } else if(jmod==9) {
    wtx <- c(0, 0, 0, 1, 1, 0, 0, 0, 0, 0.5, 0.5, 0)
    pkt <- 5 / 12
  } else if(jmod==10) {
    wtx <- c(0, 0, 0, 1, 1, 0, 0, 0, 0.5, 0.5, 0, 0)
    pkt <- 5 / 12
  } else if(jmod==11) {
    wtx <- c(0.5, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0.5)
    pkt <- 5 / 12
  } else if(jmod==12) {
    wtx <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0.5, 0.5)
    pkt <- 5 / 12
  } else if(jmod==13) {
    wtx <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0.5, 0.5, 0)
    pkt <- 5 / 12
  } else if(jmod==14) { 
    wtx <- c(0, 0, 0, 0, 1, 0, 0, 0, 0.5, 0.5, 0, 0)
    pkt <- 5 / 12
  } else { 
      stop("jmod value of ", jmod, " is invalid, must be an integer 1 to 14.")
  }
  
  rho <- exp(-phi)
  z0xx <- rho^(txx)
  r12 <- rho^(-del * c(1:12))
  con <- wtx[1] * (r12[1] - 1)
  for(k in 2:12) {
    con <- con + wtx[k] * (r12[k] - r12[k - 1])
  }
  con <- rho / (1 - rho) * con
  z0xx <- z0xx * con
  zmat <- matrix(nrow=length(txx), ncol=12)
  pckm <- matrix(nrow=length(txx), ncol=12)
  ntot <- length(txx)
  pckm[,1] <- c(txx <= del)
  for (k in 2:12) {
    pckm[,k] <- c(txx > (k - 1) * del & txx <= k * del)
  }
  zmat[,1] <- rho^txx * (rho^(-replace(txx, txx > 1/12, 1/12)) - 1)
  for (k in 2:12) {
    ztmp <- rep(0, ntot)
    for (j in 1:(k - 1)) {
      ztmp[pckm[,j]] <- 0
    }
    ztmp[pckm[,k]] <- 1 - rho^(txx[pckm[,k]] - (k - 1) * del)
    if(k < 12) {
      for (j in (k + 1):12) {
        ztmp[pckm[,j]] <- rho^(txx[pckm[,j]] - del * k) - 
          rho^(txx[pckm[,j]] - (k - 1) * del)
      }
    }
    zmat[,k] <- ztmp
  }
  sst <- z0xx
  for (k in 1:12) { 
    sst <- sst + wtx[k] * zmat[,k]
  }
  sst <- sst/phi
  medxx <- (max(sst) + min(sst)) / 2
  rngxx <- 2 * (max(sst) - medxx)
  sst <- (sst - medxx) / rngxx
  if(cmaxt <= pkt) {
    txx2 <- (txx + 1 - pkt + cmaxt)
  } else if (cmaxt > pkt)  {
    txx2 <- txx - pkt + cmaxt
  }
  txx2[txx2 > 1] <- txx2[txx2 > 1] - 1
  otmp <- order(txx2)
  sst <- sst[otmp]
  sst
}