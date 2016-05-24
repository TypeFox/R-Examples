"calcD" <- function (kincoeffs, times, cstart, rsmatrix, smatrix) {
## contributed by David Nicolaides <dnicolaides@accelrys.com>
## in modified form, along with other options for 2nd order modeling
  ns <- length(cstart)
  np <- length(kincoeffs)
  
  ## we need to pass our kinetic model the stoichiometry matrices,
  ## as well as the kinetic coefficients
  parms <- c(rsmatrix, smatrix, kincoeffs, ns, np)
  
  Model <- function(times, conc, parms) { 
    ## get s-matrices as matrices first
    rsm <- rsmatrix
    dim(rsm) <- c(ns, np)
    rsm <- t(rsm)
    sm <- smatrix
    dim(sm) <- c(ns, np)
    sm <- t(sm)
    ## now calculate the nu vectors
    nu <- rep(0,np)
    ## Note: R correctly handles 0^0 (= 1)
    for (j in 1:np) {
      nu[j] <- kincoeffs[j]*prod(conc^rsm[j,])
    }
    
    ## then calculate the sums over the nu variables to get the derivatives
    result <- rep(0,ns)
    for (j in 1:np) {
      result <- result + sm[j,]*nu[j] 
    }
    list(result) 
  }
  
  c.temp <- lsoda(y=cstart,times=times,func=Model,parms=parms)
  c.temp <- c.temp[,-1]
  colnames(c.temp) <- vector()

  c.temp
}
