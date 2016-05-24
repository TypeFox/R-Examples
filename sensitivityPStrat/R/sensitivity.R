## Function to calculate w given by equation
##
##           1
## ----------------------
##   - beta y - alpha
## %e                 + 1
.calc.w <- function(alpha, beta.y)
  1L/(1L + exp(-alpha - beta.y))

## Need to solve equaition 1 for alpha
##
##  n
## ====            dF
## \                 i
##  >    ----------------------- - C                                         (1)
## /       - beta y  - alpha
## ====            i
## i = 1 %e                  + 1
##
## to do so find the minimum of the integral of the equation 1 shown in
## equation 2
##
##  n
## ====            beta y  + alpha
## \                     i
##  >    dF  log(%e                + 1) - alpha C                            (2)
## /       i
## ====
## i = 1
##
.alpha.est <- function(alpha, beta.y, dF, C)
  sum(log(1L + exp(alpha + beta.y)) * dF) - alpha*C

.calc.alphahat <- function(beta.y, dF, C, interval) {
  alphahat <- optimize(f=.alpha.est, interval=interval, beta.y=beta.y,
                       dF=dF, C=C)$minimum
  
  if(alphahat > max(interval) - 10 || alphahat < min(interval) + 10) {
    warning("optimize overflow alphahat value invalid, alphahat = ", alphahat)
  }
  
  alphahat
}

.calc.ecdf <- function(x) {
  n <- length(x)
  
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")

  vals <- sort(unique(x))

  index <- match(x, vals)
  val.count = tabulate(index)
  return(list(F=cumsum(val.count)/n, vals=vals, index=index))
}

.foldUpperTri <- function(x) {
  xrow <- row(x)
  xcol <- col(x)
  xnrow <- nrow(x)
  upper <- xrow < xcol

  x[(xrow[upper] - 1)*xnrow + xcol[upper]] <- x[upper]
  return(x)
}

".sumCrossUpperTri<-" <- function(x, na.rm=TRUE, value) {
  xrow <- row(x)
  xcol <- col(x)
  indx <- xrow <= xcol
  
  x[indx] <- rowSums(value[xcol[indx],]*value[xrow[indx],], na.rm=na.rm)
  x
}

.makeBootstrapEvntIndx <- function(s, indx.seq, N) {
  numEvents <- 0L
  index <- integer(0L)
  storage.mode(indx.seq) <- "integer"
  storage.mode(N) <- "integer"
  
  repeat {
    subIndx <- sample(indx.seq, 1000L, replace=TRUE)
    numSubEvents <- sum(s[subIndx])
    newNumEvents <- numEvents + numSubEvents

    if(newNumEvents > N)
      subIndx <- subIndx[cumsum(s[subIndx]) + numEvents <= N]

    index <- c(index, subIndx)

    if(newNumEvents >= N)
      break

    numEvents <- newNumEvents
  }

  return(index)
}

.makeBootstrapLenIndx <- function(s, indx.seq, N)
  sample(indx.seq, N, replace=TRUE)

.calcPiPhiPsi <- function(Pi=Pi, phi=phi, psi=psi, p0=p0, p1=p1) {
  if(!missing(psi) && !is.null(psi)) {
    return(list(sens.var = "psi",
                Pi = Pi <- ifelse(psi == Inf, p1, 
                  ifelse(abs(psi) < sqrt(.Machine$double.eps), p0*p1,
                         -(sqrt((p1^2-2*p0*p1+p0^2)*exp(2*psi)+p1^2
                                +exp(psi)
                                *(-2*p1^2+2*p1-2*p0^2+2*p0)
                                +(2*p0-2)*p1+p0^2-2*p0+1)
                           +p1+exp(psi)*(-p1-p0)+p0-1)
                         /(2*exp(psi)-2))),
                psi=psi,
                phi = ifelse(psi == Inf, 1, Pi/p1)))
  } else if(!missing(phi) && !is.null(phi)) {
    new.phi <- ifelse(phi == 1, NA, phi)
    return(list(sens.var="phi",
                Pi = p1*phi,
                psi = ifelse(phi == 1,
                  Inf,
                  log((p1 * new.phi^2 + (1 - p0 - p1)*new.phi)/
                      (p1 * new.phi^2 - (p1 + p0)* new.phi + p0))),
                phi=phi))
  } else {
    return(list(sens.var = "Pi",
                Pi=Pi,
                psi = log(Pi * (1 - p1 - p0 + Pi)/(p1 - Pi)/(p0 - Pi)),
                phi = Pi/p1))
  }
}
