###### start values function ##############################################
skewhypFitStart <- function (x, breaks = NULL, startValues = "LA",
                             paramStart = NULL)   {

    if (startValues == "US") {
    #user supplied starting values
        x <- as.numeric(na.omit(x))
        svName <- "User Specified"
        if (is.null(paramStart))  stop("paramStart must be specified")
        if (!is.null(paramStart)) {
            #check parameters
            parResult <- skewhypCheckPars(paramStart)
            case <- parResult$case
            errMessage <- parResult$errMessage
            if(case == "error") stop(errMessage)
            paramStart <- c(beta = paramStart[1], delta = log(paramStart[2]),
                            mu = paramStart[3], nu = log(paramStart[4]))
        }
        names(paramStart) <- c("mu", "log(delta)","beta" , "log(nu)")
        histData <- hist(x, plot = FALSE, right = FALSE)
        if (is.null(breaks)) {
            breaks <- histData$breaks
        }
        midpoints <- histData$mids
        empDens <- ifelse(!is.finite(log(histData$density)), NA,
                          histData$density)
        maxIndex <- order(empDens, na.last = FALSE)[length(empDens)]

        out <- list(paramStart = paramStart, breaks = breaks,
                    midpoints = midpoints, empDens = empDens, svName = svName)
    }

    if(startValues == "LA"){
    #combination of moments and linear approx to log density
        x <- as.numeric(na.omit(x))
        svName <- "Linear Approximation"
        start <- skewhypFitStartLA(x, breaks = breaks)
        param <- c(start$param[1], log(start$param[2]), start$param[3],
                   log(start$param[4]))
        names(param) <- c("mu", "log(delta)","beta" , "log(nu)")
        out <- list(paramStart = param, breaks = start$breaks,
                    midpoints = start$mids, empDens = start$dens,
                    svName = svName)
    }
    return(out)
}

###### Linear Approximation function #######################################
skewhypFitStartLA <- function(x, breaks = NULL){

    #fit using combination of method of moments and a linear approxination
    #to the log density in the tails

    #histogram
    ifelse(is.null(breaks), breaks <- 30, breaks <- breaks)
    histInfo <- hist(x, plot = FALSE, breaks = breaks)
    mids <- as.matrix(histInfo$mids)
    #density
    density <- density(x)
    #find the closest density estimate to the midpoint
    func <- function(x) which.min(abs(density$x - x))
    estx <- apply(mids, 1, func)
    esty <- density$y[estx]
    #lower and upper tail
    quantiles <- quantile(x, probs = c(0.1,0.9))
    lower <-  which(mids <= quantiles[1])
    upper <-  which(mids >= quantiles[2])
    xLower <- mids[lower]
    xUpper <- mids[upper]
    yLower <- esty[lower]
    yUpper <- esty[upper]
    #stop if either tail has less than two points
    if( length(xLower) < 4 | length(xUpper) < 4 )
        stop("not enough points to model tail behaviour")

    #determine skewness
    skew <- skewness(x)

    ##SYMMETRY ##
    if( abs(skew) < 0.1 ) {
        beta <- 0
        mu <-  mean(x)
        #fit models
        upperModel <- lm( log(yUpper) ~ log(abs(xUpper)) )
        lowerModel <- lm( log(yLower) ~ log(abs(xLower)) )
        #estimates
        nuLow <- -2*(upperModel$coeff[2] + 1)
        nuUpp <- -2*(lowerModel$coeff[2] + 1)
        #average to get nu
        nu <- mean(c(nuLow,nuUpp))
        #nu must be >4
        if(nu <= 4) nu <- 4.1
        #find delta
        delta <- sqrt( var(x)*(nu - 2) )
    }

    ## POSITIVE SKEW ##
    if( skew >= 0.1 ){#heavy tail is upper tail
        #fit models
        lowerModel <- lm( log(yLower) ~ log(abs(xLower)) + abs(xLower) )
        upperModel <- lm( log(yUpper) ~ log(abs(xUpper)) )
        #estimates
        estLow <- lowerModel$coeff[c(2,3)]
        nuLow <- -2*(estLow[1] + 1)
        betaLow <- (estLow[2]/ - 2)
        estUpp <- upperModel$coeff[2]
        nuUpp <- -2*(estUpp + 1)
        #nu must be > 0
        ifelse(nuLow > 0, nu <- nuLow, nu <- nuUpp)
        beta <- betaLow
        #nu must be > 4
        if(nu <= 4) nu <- 4.01
        #solve for delta
        var <- var(x)
        sol1 <- - 1/4*(nu - 4 - (nu^2 - 8*nu + 16 + 8*beta^2*var*nu -
                                 32*beta^2*var)^(1/2))/beta^2*(nu - 2)
        sol2 <- -1/4*(nu - 4 + (nu^2 - 8*nu + 16 + 8*beta^2*var*nu -
                                32*beta^2*var)^(1/2))/beta^2*(nu - 2)
        #delta^2 must be positive
        ifelse(sol1 > 0, delta2 <- sol1, delta2 <- sol2)
        delta <- sqrt(delta2)
        #solve for mu
        mu <-  mean(x) - beta*delta2/(nu - 2)
    }

    ## NEGATIVE SKEW ##
    if( skew <= -0.1){ #heavy tail is lower tail
        #fit models
        lowerModel <- lm( log(yLower) ~ log(abs(xLower)) )
        upperModel <- lm( log(yUpper) ~ log(abs(xUpper)) + abs(xUpper) )
        #estimates
        estLow <- lowerModel$coeff[2]
        nuLow <- -2*(estLow + 1)
        estUpp <- upperModel$coeff[c(2,3)]
        nuUpp <- -2*(estUpp[1] + 1)
        betaUpp <- -(estUpp[2] / -2)
        #nu must be > 0
        ifelse(nuLow > 0, nu <- nuLow, nu <- nuUpp)
        beta <- betaUpp
        #nu must be > 4
        if(nu <= 4) nu <- 4.1
        #solve for delta
        var <- var(x)
        sol1 <- -1/4*(nu - 4 - (nu^2 - 8*nu + 16 + 8*beta^2*var*nu -
                                32*beta^2*var)^(1/2))/beta^2*(nu - 2)
        sol2 <- -1/4*(nu - 4 + (nu^2 - 8*nu + 16 + 8*beta^2*var*nu -
                              32*beta^2*var)^(1/2))/beta^2*(nu - 2)
        #delta^2 must be positive
        ifelse(sol1 > 0, delta2 <- sol1, delta2 <- sol2)
        delta <- sqrt(delta2)
        #solve for mu
        mu <- mean(x) - beta*delta2/(nu-2)
    }

    #return results
    param <- c(mu,delta,beta,nu)
    names(param) <- c("mu","delta","beta","nu")
    out <- list(param = param, breaks = histInfo$breaks,
                mids = as.vector(mids), dens = esty)
    return(out)
}
