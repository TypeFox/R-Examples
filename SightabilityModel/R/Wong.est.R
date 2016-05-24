Wong.est <-
function(total, srates, nh, Nh, stratum, subunit, covars, beta, varbeta, smat=NULL){

  #  Arguments (the first four should have one observation for each individual in the obs. dataset
  # total = total number of animals in the observed group
  # srates = sampling rate for the stratum (for each observed animal)
  # nh = number of sampled units in each stratum
  # Nh = number of population units in each stratum
  # stratum = stratum identifier  
  # Subunit = subunit identifier 
  # covars = matrix of covariates used in the sightability model
  # beta = vector of parameter estimates from logistic model
  # varbeta = var/cov matrix for above parameter estimates
  # smat = variance covariance matrix of estimated inflation factors (theta)
  
  # Check on length of vector arguments
    n <- length(total)
    if(length(srates) != n) {stop("Srates vector needs to be same dimension as total vector")}
    if(length(stratum) != n) {stop("Stratum vector needs to be same dimension as total vector")}
    if(length(srates) != n) {stop("Subunit vector needs to be same dimension as total vector")}

  # Form xmatrix
    xdata <- as.matrix(cbind(rep(1, n), covars))
    ncovars <- dim(xdata)
    if(ncovars[1] != n) {stop("Covariate matrix needs to be same length as total vector")}  
      
  # Make sure beta is a matrix of dim ncovar x 1 
    beta <- matrix(beta, length(beta), 1)
    if (length(beta) != ncovars[2]){stop("Covars matrix must have same number of columns as the Beta vector")} 
 
  #  Estimate theta for each animal = 1/prob(detection)
    theta <- 1+exp(-xdata%*%beta-diag(xdata%*%varbeta%*%t(xdata))/2)  

  # Xie = prob. detection = 1/ theta  
    xie <- 1/theta
  
  # Point estimate
    tau.hat <- sum(total*theta/srates) 
 
  #------------------------- Now, calculate the variance components in steps  ----------------------------------
  #  All 3 variance components require calculation of Cov(theta^[ij], theta^[i'j']) (see Wong 1996)
  #  So, start by calculating variance covariance matrix of theta values = smat (below)
    nxdata <- nrow(xdata)
    sm3 <- matrix(0, nxdata, nxdata) # holder for some of the terms in the expression for smat
   
   if(is.null(smat)){ #Variance/covariance matrix of smat not supplied by user
                    # Use asymptotic estimate based on log-normal assumption (of SS and Wong)  
  # Do as much of the matrix multiplication outside of loop as possible
    xb <- xdata%*%beta  # X'beta
    xbb <- kronecker(xb, t(xb), FUN = "+")  # x1+x2
    xvarbeta <- xdata%*%varbeta%*%t(xdata)  # X Sig X
    for(i in 1:nxdata){
      for(j in i:nxdata){
       xtemp1 <- as.vector(xdata[i, ], mode = "numeric")
       xtemp2 <- as.vector(xdata[j, ], mode = "numeric")
       xtot <- t(xtemp1+xtemp2)
       sm3[i, j] <- sm3[j, i] <- xtot%*%varbeta%*%t(xtot)/2
      }
    }
    smat <- exp(-xbb-sm3)*(exp(xvarbeta)-1)
  }
  
  #------------------------ End Calculation of Cov matrix ------------------------------------------------

  #---------------------- Calculate Var(T|D) = "Sampling variance" ---------------------------------------

  #  Start by calculating eq. 2.2.4 in Wong (1996)

  #  Calculate MKs and pks (MKs = corrected totals by subunit, pks = sampling rate for the subunit)
    MKs <- rowsum(total*theta, subunit)
    pks <- tapply(srates,subunit, unique)
    nk <- length(pks)
 
  # Turn MKs and pks into vectors 
    MKs <- matrix(MKs, length(MKs), 1)
    pks <- matrix(pks, length(pks), 1)

  # Calculate second term in Eq 2.2.4 of Wong
  # kdata = (stratum id, subunit ID)
    temp <- tapply(stratum, subunit, unique) # just pulls off stratum ids
    kdata <- data.frame(cbind(temp, as.numeric(names(temp))))
    names(kdata) <- c("stratum", "k")
  
  # Re-number stratum so that go from 1 to h (if not already coded this way)
    kdata[, 1] <- as.numeric(factor(kdata[, 1]))  

  # pkkprime for 2 observations from same stratum,h = (nh/Nh)*(nh-1)/(Nh-1)
    pkkprime <- as.vector(nh*(nh-1)/(Nh*(Nh-1)))

  # Update to fix bug noted by Cliff Rice (when applied to a single sampling unit)
   if(all(nh == Nh) != TRUE){

  # Now, calculate the sum [(pkk'-pk*pk')/pkk'*pk*pk']*Mk*MK'
    sterm <- 0

  # Update 2-7-2012:  fix bug if there is only 1 sample plot with observed animals...then,
  #  no covariance terms below
    if(nrow(kdata) > 1){

       for(i in 1:(nk-1)){
         for(j in (i+1):nk){

  #  if the two observations are in different strata, then pkk'=pk*pk' => adds 0 to sterm
             if(kdata[i, 1] == kdata[j, 1]){ # same stratum
              sterm <- sterm+
              (pkkprime[kdata[i, 1]]-pks[i]*pks[j])*MKs[i]*MKs[j]/(pkkprime[kdata[i, 1]]*pks[i]*pks[j])}
             }
         }
       } 
       Eq2.2.4 <- sum(((1-pks)/pks^2)*MKs^2)+2*sterm
  
    }else{Eq2.2.4 <- 0}

  # Now, calculate the other terms that need to be subtracted:
    Eq2.2.7t3 <- sum(((1-srates)/srates^2)*total^2*(theta^2-theta))
  
  #  Next 2 terms are easier to calculate w/ loops (same w/ similar components of model variance)
    Eq2.2.7t4 <- 0
    Eq2.2.7t5 <- 0
    Eq2.2.12t2 <- 0
    Eq2.2.12t3 <- 0
    for(i in 1:nxdata){
      for(j in 1:nxdata){
        if(subunit[i] == subunit[j]){
          temp <-(1/srates[i]^2)*total[i]*total[j]*smat[i, j]
          Eq2.2.7t4 <- Eq2.2.7t4+(1-srates[i])*temp*I(i != j)
          Eq2.2.12t2 <- Eq2.2.12t2+temp*I(i != j)
         }
  # NOte last term of Eq 2.2.7 is 0 if the two observations are from the same subunit (i=i') or,
  #    if they come from different strata (because then pii'=pi*pi')  
  # But, last term of EQ2.2.12 is not 0 if come from different strata  
        if(subunit[i] != subunit[j]){
          if(stratum[i] == stratum[j]){
             Eq2.2.7t5<-Eq2.2.7t5+((pkkprime[stratum[i]]-srates[i]*srates[j])/(pkkprime[stratum[i]]*srates[i]*srates[j]))*total[i]*total[j]*smat[i, j]
          }
          Eq2.2.12t3 <- Eq2.2.12t3+(1/(srates[i]*srates[j]))*total[i]*total[j]*smat[i, j]
        }
      }   
    } 
     
    VarSamp <- Eq2.2.4- Eq2.2.7t3 - Eq2.2.7t4 -  Eq2.2.7t5
   
  # Now, E(var(T|D)) = sightability variability
    VarSight <- sum((total^2/srates^2)*(theta^2-theta-diag(smat)))
  
  # Now, model component
    m1 <- sum((total^2/srates^2)*diag(smat))
    Varmod <- m1+Eq2.2.12t2+Eq2.2.12t3

    Vartot <- VarSamp+VarSight+Varmod

    out <- c(tau.hat, Vartot, VarSamp, VarSight, Varmod)
    names(out) <- c("tau.hat", "VarTot", "VarSamp", "VarSight", "VarMod")
    out2 <- NULL
    out2$est <- out
    out2$var.method = "WONG"
    out2
}
