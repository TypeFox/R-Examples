SS.est <-
function(total, srates, nh, Nh, stratum, subunit, covars, beta, varbeta, smat = NULL){

  # Arguments (the first four should have one observation for each individual in the obs. dataset
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
    beta < -matrix(beta, length(beta), 1)
    if (length(beta) != ncovars[2]){stop("Covars matrix must have same number of columns as the Beta vector")} 
 
  #  Estimate theta for each animal = 1/prob(detection)
    theta <- 1+exp(-xdata%*%beta-diag(xdata%*%varbeta%*%t(xdata))/2)
  
  # Xie = pie = prob. detection = 1/ theta 
  # Note: an alternative (the MLE) is inv.logit(exp(xB) instead of 1/theta), but in simulatons, 1/theta seems to work better
    xie <- 1/theta
   
  # Point estimate
    tau.hat <- sum(total*theta/srates)
 
  # Now, calculate the variance in steps 
  # First, calculate Var(T|D) = "Sampling variance" (v1 of Thopmson and Seber 1996)

  #  Calculate MKs and pks (MKs = corrected totals by subunit, pks = sampling rate for the subunit)
    MKs <- rowsum(total*theta, subunit)
    pks <- tapply(srates,subunit, unique)
    nk <- length(pks)
 
  # Turn MKs and pks into vectors 
    MKs <- matrix(MKs, length(MKs), 1)
    pks <- matrix(pks, length(pks), 1)

  # Calculate second term in varT.D
  # kdata = (subunit, stratum id)
    temp <- tapply(stratum, subunit, unique) # just pulls off stratum ids
    kdata <- data.frame(cbind(temp, as.numeric(names(temp))))
    names(kdata) <- c("stratum", "k")
  
  # Renumber stratum so that go from 1 to nh
    kdata[ , 1] <- as.numeric(factor(kdata[ , 1]))  

  # pkkprime for stratum 1, 2, 3,...  = (nh/Nh)*(nh-1)/(Nh-1)
    pkkprime <- as.vector(nh*(nh-1)/(Nh*(Nh-1)))
  
  # Update to fix bug noted by Cliff Rice (when applied to a single sampling unit)
    if(all(nh == Nh) != TRUE){
  
    # Counter for variance component
    sterm <- 0

  # Update 2-7-2012:  fix bug if there is only 1 sample plot with observed animals...then,
  #  no covariance terms below
    if(nrow(kdata) > 1){

      for(i in 1:(nk-1)){
        for(j in (i+1):nk){
        #  if the two observations are in different strata, then pkk'=pk*pk' => adds 0 to sterm
          if(kdata[i, 1] == kdata[j, 1]){ # same stratum
            sterm<-sterm+
            (pkkprime[kdata[i, 1]]-pks[i]*pks[j])*MKs[i]*MKs[j]/(pkkprime[kdata[i, 1]]*pks[i]*pks[j])}
          }
        }
       } 
    varT.D <- sum(((1-pks)/pks^2)*MKs^2)+2*sterm
    }else{varT.D <- 0}

  # Now, E(var(T|D)) = "sightability variability" (v2 of Thompson and Seber 1996)
  # Note:  original paper by Steinhorst & Samuel had srates^2, but paper in 1992 has srates
  #          and Thompson suggests using srates.
  #   EvarT.D<-sum(((1-xie)/xie^2)*total^2/srates^2) = vd (from Steinhorst and Samuel 1989)
    EvarT.D <- sum(((1-xie)/xie^2)*total^2/srates)  #this version is better, but earlier version to compare w/ older versions

  #  Now, calculate uncertainty due to estimating thetas 
    uxdata2 <- as.matrix(unique(xdata))
    nxdata <- nrow(uxdata2)
    sm3 <- matrix(0, nxdata, nxdata)
   
  # Do as much of the matrix multiplication outside of loop as possible
    xb <- uxdata2%*%beta
    xbb <- kronecker(xb, t(xb), FUN = "+")
    xvarbeta <- uxdata2%*%varbeta%*%t(uxdata2)  
    for(i in 1:nxdata){
      for(j in i:nxdata){
        xtemp1 <- as.vector(uxdata2[i,], mode = "numeric")
        xtemp2 <- as.vector(uxdata2[j,], mode = "numeric")
        xtot<-t(xtemp1+xtemp2)
        sm3[i, j] <- sm3[j, i] <- xtot%*%varbeta%*%t(xtot)/2
      }
    }
    smat <- exp(-xbb-sm3)*(exp(xvarbeta)-1)
  
  # Determine ajs in SS's estimator of Vm
    index2 <- as.numeric(factor(apply(xdata, 1, paste, collapse = ":"), 
                  levels=unique(apply(xdata, 1, paste, collapse = ":"))))  
    aj <- tapply(total/srates, index2, sum)
  
    Var.mod <- t(aj)%*%smat%*%aj 
    vartot <- varT.D+EvarT.D+Var.mod
    out <- c(tau.hat, vartot, varT.D, EvarT.D, Var.mod)
    names(out) <- c("tau.hat", "VarTot", "VarSamp", "VarSight","VarMod")
    out2 <- NULL
    out2$est <- out
    out2$var.method = "SS"
    out2
}
