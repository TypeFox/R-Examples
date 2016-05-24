## transfer moments about different locations for any distributions
momChangeAbout <- function(order = "all", oldMom, oldAbout, newAbout) {
     if (!is.vector(oldMom)){
        stop("A vector of moments must be supplied")
     }
     if (order == "all") {
        ## Compute moment of up to length(oldMom) about location new
        mom <- rep(NA,length(oldMom))
        oldMoment <- c(1,oldMom)
        for (i in 1:length(oldMom)) {
          oldMom <- oldMoment[1:(i+1)]
          binomCoeff <- choose(i, 0:i)
          diffPower <- (oldAbout - newAbout)^(i:0)
          mom[i] <- sum(binomCoeff*diffPower*oldMom)
        }
     } else {
       ## Check order is within in the right range
       if (length(oldMom) < order) {
          stop("The length of of the vector oldMom must not be less than the
                value of order")
       }
       if (!is.wholenumber(order)){
          stop("Order must be a whole number")
       }
       if ((order < 0)) {
          stop("Order must be positive")
       }
       ## Compute moment of a specific order about location new
       oldMom <- c(1,oldMom)
       oldMom <- oldMom[1:(order+1)]
       binomCoeff <- choose(order, 0:order)
       diffPower <- (oldAbout - newAbout)^(order:0)
       mom <- sum(binomCoeff*diffPower*oldMom)
     }

     ## Return moment
     return(mom)
}

