
################################################################################
#######################  AUTOCORRELATION AND PARTIAL ###########################
################################################################################

z.normalize = function(x) {
    (x - mean(x)) / sd(x)
}

range.normalize = function(x) {
    minim <- min(x)
    maxim <- max(x)
    (x -minim) / (maxim - minim)    
}

.common.ts.sanity.check <- function(x) {
    if (missing(x)) {
        stop("At least one series is missing!")
    }
    if (any(is.na(x))) {
        stop("NA in the series")
    }
    if (!is.numeric(x)) {
        stop("Series must be numeric")
    }
    #check length
    if (length(x) < 2) {
        stop("Incorrect length of the series")
    }
    if (!is.null(dim(x))) {
        stop("Incorrect dimension of the series, please input univarate series")
    }
}
.ts.freq.check <- function(x, y) {
    if (is.ts(x) && is.ts(y)) { #check their frequencies
        cbind(x,y)
    }
}
.ts.sanity.check <- function(x,y) {
    .common.ts.sanity.check(x)
    .common.ts.sanity.check(y)  
    .ts.freq.check(x,y)
}

#check if series have equal length, a requisite of some functions
.check.equal.length.ts <- function(x,y) {
    if (length(x) != length(y)) {
        stop("Time series must have the same length")
    }	
    
}

#weighted distance of acf and pacf coefficients
.internal.autocorr.dist <- function(rhox, rhoy, p=NULL, omega=NULL) {
    if ( length(rhox) != length(rhoy) ) {#check compatible coefficient vectors
        stop("The amount of autocorrelation coefficients must be the same, maybe lag.max greater than the length of one of the series")
    }
    if (is.null(omega)) { #if there is no weighting matrix
        if (!is.null(p)) { #check if there is gemoetrical decay parameter
            omega <- diag(p*(1-p)**(1:length(rhox)))
        }
        else { #no weightinh matrix and no geomtrical decay parameter, use identity matrix
            omega <- diag(length(rhox))
        }
    }
    t(rhox - rhoy) %*% omega %*% (rhox - rhoy) #weighted euclidean distance
}


diss.ACF <- function(x, y ,  p=NULL,  omega=NULL, lag.max=50) {
    .ts.sanity.check(x, y)
    rhox <- acf(x, lag.max=lag.max, plot=FALSE)$acf[-1]
    rhoy <- acf(y, lag.max=lag.max, plot=FALSE)$acf[-1]
    .internal.autocorr.dist( rhox, rhoy, p, omega)
}



diss.PACF <- function(x, y, p = NULL, omega=NULL, lag.max=50) {
    .ts.sanity.check(x, y)
    rhox <- as.vector(pacf(x, lag.max=lag.max, plot=FALSE)$acf)
    rhoy <- as.vector(pacf(y, lag.max=lag.max, plot=FALSE)$acf)
    .internal.autocorr.dist( rhox, rhoy, p, omega)
}





#######################################################
##########  distance Piccolo  #########################
#######################################################


#try to find the ar coefficients of a AR series, if no order found, try forcing
find_ar_model_force = function(x, permissive) {
    arx <- ar(x) #first, try to find automatically
    if (arx$order < 1) { #if no order found, try forcing order 1
        if (permissive) {
            arx <- ar(x, aic=FALSE, order.max = 1)
        }
        if (arx$order < 1) stop("Could not find a valid AR order for the series")
    }
    arx
}


#######################################

diss.AR.PIC <- function(x, y, order.x=NULL, order.y=NULL, permissive=TRUE) {
    .ts.sanity.check(x, y)
    #if order NULL use ar AIC, else use arima fitting of the given order
    PIx <- NULL
    if (is.null(order.x)) { #no ARIMA order, use AR
        PIx <- find_ar_model_force(x, permissive)$ar
    } else { 
        if ((order.x[1]) < 1) stop("The arima order must have AR coefficients, they are used for the distance")
        arim <- arima(x, order.x)
        PIx <- arim$coef[1:order.x[1]] #get the AR coeff off the arima model
    }
    
    PIy <- NULL
    if (is.null(order.y)) { #no ARIMA order, use AR
        PIy <- find_ar_model_force(y, permissive)$ar
    } else {
        if ((order.y[1]) < 1) stop("The arima order must have AR coefficients, they are used for the distance")
        arim <- arima(y, order.y)
        PIy <- arim$coef[1:order.y[1]] #get the AR coeff off the arima model
    }
    
    k <- max(c(length(PIx), length(PIy))) #get the maximun order
    if (k < 1) {
       stop("Could not find any AR coefficients")
    }
    
    PRIMAx <- rep(0,k) #fill with zeroes to the greatest AR order between series x and y (k)
    if (length(PIx) > 0) {
        PRIMAx[1:length(PIx)] <- PIx
    }
    
    PRIMAy <- rep(0,k)
    if (length(PIy) > 0) {
        PRIMAy[1:length(PIy)] <- PIy
    }
    
    as.numeric( dist(rbind(PRIMAx,PRIMAy)) ) #compute the euclidean distance between the zero padded AR coefficients
}


##################################################
########### distance Maharaj  ####################
##################################################

#regression model
maharajahextended <- function( x, k ) {
    
    X <- x[-(1:k)]
    
    TT <- length(x)
    
    Wx <- matrix(ncol=k, nrow=(TT-k))
    
    for (i in 1:(TT -k) ) {
        Wx[i,] <- x[(k +i - 1):(i)]
    }
    result <- list()
    result$X <- X
    result$Wx <- Wx
    result
}

#extended distance see reference article in the documentation
distance.MAH.EXT <- function( x, y, k) {
    MX <- maharajahextended(x, k)
    MY <- maharajahextended(y, k)
    w <- dim(MX$Wx)[2]
    h <- dim( MX$Wx)[1]
    bigW <- matrix(0, nrow=2*h, ncol=2*w)
    for ( j in 1:w) {
        for (i in 1:h) {
            bigW[i,j] <- MX$Wx[i,j]  
        }
    }
    for ( j in 1:w) {
        for (i in 1:h) {
            bigW[h+i,w+j] <- MY$Wx[i,j]  
        }
    }
    Epsil <- cov(cbind(x,y))
    bigW
    Epsil
    Iden <- diag(1, length(x) - k)
    Iden
    V <- kronecker(Epsil, Iden  )
    IV <- solve(V)
    tryCatch( { #these operations can go wrong, it they fail, try a smaller order
        PI <- solve(t(bigW) %*% IV %*% bigW) %*% t(bigW) %*% IV %*% c(MX$X, MY$X)
        R <- cbind( diag(1, k), diag(-1,k) )
        result <- list()
        result$statistic <- t(R %*% PI) %*% solve( R %*% (t(bigW) %*% IV %*% bigW ) %*% t(R)) %*% ( R%*%PI)
        result$p_value <- pchisq(result$statistic, k, lower.tail=F)
        result
    }, error = function(e) { if (k>1) {
        distance.MAH.EXT(x,y,k-1)
    }
      else { stop("Could not find valid AR order")}
    })
}



distance.MAH.SIMP = function( x, y, PIx, PIy, var.pred.x, var.pred.y, k, permissive=TRUE ) {
    PRIMAx <- rep(0,k)   #fill with zeroes
    if (length(PIx) > 0) {
        PRIMAx[1:length(PIx)] <- PIx
    }
    PRIMAy <- rep(0,k)
    if (length(PIy) > 0) {
        PRIMAy[1:length(PIy)] <- PIy
    }
    
    covx <- acf(x, lag.max=k-1, type="covariance", plot=FALSE)$acf
    covy <- acf(y, lag.max=k-1,type="covariance", plot=FALSE)$acf
    Rx <- matrix(nrow=length(covx),ncol=length(covx))
    Ry <- Rx
    for (i in 1:length(covx) ) {
        indices <- (c(i:1, 2:(length(covx)-i+1)))
        Rx[i,] <- covx[ indices[1:length(covx)]]
        Ry[i,] <- covy[ indices[1:length(covx)]]
    }
    V <- (solve(solve(Rx)*(var.pred.x) + solve(Ry)*(var.pred.y)))
    dif <- (PRIMAx - PRIMAy)
    
    D <- sqrt(length(x)) * ( dif %*% V %*% dif)
    
    list(statistic=D, p_value=pchisq(D, k, lower.tail=F))
}


diss.AR.MAH = function( x, y, dependence = FALSE, permissive = TRUE) {
    .ts.sanity.check(x, y)
    arx <- find_ar_model_force(x, permissive)
    ary <- find_ar_model_force(y, permissive)
    k <- max(c(length(arx$ar), length(ary$ar)))
    if (k < 1) {
        if (k < 1) stop("Could not find a valid AR order for the series")
    }

    if (dependence) {
        .check.equal.length.ts(x,y)
        distance.MAH.EXT(x, y, k)
    }
    else {
        distance.MAH.SIMP(x, y, arx$ar, ary$ar, arx$var.pred, ary$var.pred, k, permissive)
    }
}

######################################################
########  PERIODOGRAM BASED DISTANCES  ###############
######################################################

diss.PER <- function(x,y, logarithm=FALSE, normalize=FALSE) {
    .ts.sanity.check(x, y)
    .check.equal.length.ts(x,y)
    Ix <- spec.pgram(x, plot=F)$spec
    Iy <- spec.pgram(y, plot=F)$spec
    if (normalize) {
        Ix <- Ix/var(x)
        Iy <- Iy/var(y)
    }
    if (logarithm) {
        Ix <- log(Ix)
        Iy <- log(Iy)
    }
    dist(rbind(Ix,Iy))/(length(Ix))
}


diss.INT.PER <- function(x,y, normalize=TRUE) {
    .ts.sanity.check(x, y)
    .check.equal.length.ts(x,y)
    Ix <- spec.pgram(x, plot=F)
    Iy <- spec.pgram(y, plot=F)
    Cx <- 1
    Cy <- 1
    if (normalize) {
        Cx <- sum(Ix$spec)
        Cy <- sum(Iy$spec)
    }
    sum ( abs(cumsum(Ix$spec)/Cx - cumsum(Iy$spec)/Cy) )
}


################################################
###  SPECTRAL DENSITY APPROXMATION DISTANCES ###
################################################


### maximum likelihood functions ###

#kernel function
funcionKh <- function( value, h) {
    value <- value/h
    dnorm( -(value**2) ) / h
}

#function to be optimized for maximum likelihood, see referenced papers
Spectral.AB <- function (  ABvec, lambda, Yks, lambdas, h) {
    acum <- 0
    a <- ABvec[1]
    b <- ABvec[2]
    
    -sum ( ( -exp(Yks - a -b*(lambdas - lambda) ) + Yks - a -b*(lambdas - lambda) ) * funcionKh( lambdas - lambda, h) )
    
}

likelihood.optim <- function(  lambda, Yks, lambdas, h) {
    
    startA = lambdas[1]
    startB = 0
    
    optim(c(startA,startB), Spectral.AB, lambda=lambda, Yks = Yks, lambdas=lambdas, method="L-BFGS-B",lower = c(min(Yks) , -101), upper= c(max(Yks),101), h =h)$par[1]
}
### end of maximul likelihood


#diveregence function
divergenceW <- function( x, alpha) {
    if ((alpha > 0) & (alpha< 1)) {
        log(alpha*x + 1 - alpha) - alpha*log(x)
    }
    else {
        stop("condition 0 < alpha < 1 not fulfilled")
    }
}

simetricDivergenceW <- function(x,alpha) {
    divergenceW(x,alpha) + divergenceW(1/x,alpha)  
}

#plot the soothed spectral density
plotsmoothspec <- function ( lambdasX, YksX, lambdasY, YksY, myf, hX, hY, n, method="Maximum Likelihood") {
    if ( n < 1) {
        n <- 500
    }
    baseX <- seq(min(lambdasX) + 0.001, max(lambdasX) - 0.001, length.out=min(500,n) )
    specX <- NULL
    if (pmatch(method , c("Maximum Likelihood", "Least Squares")) == 2) {
        specX <- (myf(YksX, lambdasX,  hX, baseX))
    } else {
        specX <- exp(myf(baseX, YksX, lambdasX, hX))
    }
    baseY <- seq(min(lambdasY)+ 0.001, max(lambdasY)- 0.001, length.out=min(500, n) )
    specY <- NULL
    if (pmatch(method , c("Maximum Likelihood", "Least Squares")) == 2) {
        specY <- (myf(YksY,lambdasY, hY, baseY))
    } else {
        specY <- exp(myf(baseY, YksY, lambdasY, hY))
    }
    plot.default(baseX, specX,type="l", col="red", ylim=c( min(specX,specY), max(specX, specY) ) ,
                 main=paste("Approx. spectral density by ", method), xlab="frequency",ylab="spectrum")
    lines(baseY, specY,type="l", col="blue")
    legend("topright", pch=16, col=c("red", "blue"), legend= c("x", "y") )
    
    
}

.vectorized.lk.optim <- Vectorize(likelihood.optim,"lambda")  #needed for likelihood.optim to accept a vector, required for integrate


#trapezoind integration, taken from another package
trapez <- function(x,y) {
    idx = 2:length(x)
    return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}

interp.SPEC.LOGLIKELIHOOD <- function(x, n) {
    pgram <- spec.pgram(x, plot=FALSE)
    Yks <- log(pgram$spec)
    lambdas <- pgram$freq
    interplambdas <- seq(min(lambdas), max(lambdas), length.out=n)
    hX <- 0.93*dpill(lambdas, Yks)
    approx( interplambdas, .vectorized.lk.optim( interplambdas, Yks, lambdas, hX ) )
}

interp.W.LK <- function(x, n) {
    interps <- interp.SPEC.LOGLIKELIHOOD(x, n)
    interps$y <- exp(interps$y)
    interps
}

integrate.divergenceW <- function(base, x, y, alpha ) {
    val <- simetricDivergenceW( x / y, alpha)
    trapez(base, val)
}

leastsquares.spec <- function( Yk, lambdas, h, lambdaeval) {
    d <- data.frame(Yk)
    d$lambdas <- lambdas
    lp <- locpol(Yk~lambdas, d, bw=h,kernel=gaussK, xeval=lambdaeval )
    ord <- order(lambdaeval) #trick to get the original order of the lambas, locpol sorts the input vector
    ord2 <- order(ord)   #second part of the trick
    lp$lpFit$Yk[ord2]
}

interp.SPEC.LS <- function(x, n)  {
    pgram <- spec.pgram(x, plot=FALSE)
    Yk <- pgram$spec
    lambdas <- pgram$freq
    h <- dpill(lambdas, Yk)
    interplambdas <- seq(min(lambdas), max(lambdas), length.out=n)
    ys <- leastsquares.spec(Yk, lambdas, h, interplambdas)
    ys[ys<0.0001] <- 0.0001 #no zeroes allowed
    list( x = interplambdas, y = ys)
}

integrate.ISD <- function(base, x, y) {
    trapez(base, (x - y)^2)
}

interp.SPEC.GLK <- function(x, n) {
    pgram <- spec.pgram(x, plot=FALSE)
    Yk <- log(pgram$spec)
    lambdas <- pgram$freq
    h <- 0.93*dpill(lambdas, Yk)
    ys <- .vectorized.lk.optim( lambdas, Yk, lambdas, h )
    list( x = lambdas, y = list(mu = ys, Z = Yk) )
}

integrate.GLK <- function( base, x, y) {
    Z <- x$Z - y$Z
    mu <- x$mu - y$mu
    sum(Z - mu - 2*log(1 + exp(Z - mu))) - sum( Z - 2*log(1 + exp(Z)))
}
#generic linear interpolation approximation of the spectral dissimilarities
#interpfun is a function to calculate the spectrum aproximation
#n the size of the grid for interpolation
#intergrationfun is calculates the sum/integration of the diferences
multidiss.interp.SPEC <- function( series, n, interpfun, integrationfun, ...) {
    l <- length(series)
    dists <- matrix(0, l, l)
    #get the interpolated values
    interps <- lapply(series, interpfun, n)
    base <- interps[[1]]$x
    ##calc the function with the interpolated values
    for (i in 1:(l-1)) {
        for (j in (i+1):l) {
            d <- integrationfun( base, interps[[i]]$y , interps[[j]]$y, ...)
            dists[i,j] <- d
            dists[j,i] <- d
        }
    }
    as.dist(dists)
}

distance.W.LK <- function(x,y, alpha, plot=FALSE, n=length(x)) {
    
    
    
    pgx <- spec.pgram(x,plot=FALSE)
    YksX <- log(pgx$spec)
    lambdasX <- pgx$freq
    pgy <- spec.pgram(y,plot=FALSE)
    YksY <- log(pgy$spec)
    lambdasY <- (pgy$freq)
    
    hX <- 0.93*dpill(lambdasX, YksX)
    hY <- 0.93*dpill(lambdasY, YksY)
    
    integrateaux <- function( lambda ) {
        xx <- exp(.vectorized.lk.optim(lambda, YksX, lambdasX, hX))
        yy <- exp(.vectorized.lk.optim(lambda, YksY, lambdasY, hY))
        simetricDivergenceW(  xx / yy, alpha)
    }
    
    a <- 0
    if (n > 0) {
        a <- multidiss.interp.SPEC(list(x,y), n, interp.W.LK, integrate.divergenceW, alpha)
    } else {
        tryCatch( {
            a <- integrate(integrateaux, min(lambdasX), max(lambdasX))$value
        }, error = function (e) {
            warning("Failed approximation with window from plug.in method, increasing window...")
            hX <- 2*hX
            hY <- 2*hY
            a <- integrate(integrateaux, min(lambdasX), max(lambdasX))$value
        })
    }
    if (plot) {
        plotsmoothspec(lambdasX, YksX, lambdasY, YksY, .vectorized.lk.optim, hX, hY, n)
    }
    a
}

distance.W.DLS <- function(x, y, alpha, plot=FALSE, n=length(x)) {
     
    YksX <- (spec.pgram(x,plot=FALSE)$spec)
    lambdasX <- spec.pgram(x,plot=FALSE)$freq
    YksY <- (spec.pgram(y,plot=FALSE)$spec)
    lambdasY <- (spec.pgram(y,plot=FALSE)$freq)
    hX <- dpill(lambdasX, YksX)
    hY <- dpill(lambdasY, YksY)
    
    integrateaux <- function( lambda ) {
        xx <- leastsquares.spec(YksX, lambdasX, hX, lambda )
        yy <- leastsquares.spec(YksY, lambdasY, hY, lambda )
        xx[xx<0.0001] <- 0.0001
        yy[yy<0.0001] <- 0.0001
        simetricDivergenceW(  xx / yy, alpha)
    }
    lambdas <- spec.pgram(x, plot=F)$freq
    
    a <- 0
    if (n > 0) {
        a <- multidiss.interp.SPEC(list(x,y), n, interp.SPEC.LS, integrate.divergenceW, alpha)
    } else {
        a <- integrate(integrateaux, min(lambdas), max(lambdas), subdivisions=100)$value
    }
    if (plot) {
        plotsmoothspec(lambdasX, YksX, lambdasY, YksY, leastsquares.spec, hX, hY, n, "Least Squares")
    }
    
    a
}


diss.SPEC.LLR <- function(x,y, alpha=0.5, method="DLS", plot=FALSE, n=length(x)) {
    .ts.sanity.check(x, y)
    .check.equal.length.ts(x,y)
    typedist <- 0
    type <-  (pmatch(method, c("DLS", "LK" )))
    if (is.na(type)) {
        stop(paste("Unknown method", method))
    } else if (type == 1) {
        typedist <- distance.W.DLS(x,y, alpha, plot, n)
    }
    else if (type == 2) {
        typedist <- distance.W.LK(x,y, alpha, plot, n)
    }
    typedist
}

multidiss.SPEC.LLR <- function(series, method="DLS", alpha=0.5, plot=FALSE, n=length(series[[1]])) {
    if ( n > 0) {
        interpfun <- NULL
        type <-  (pmatch(method, c("DLS", "LK" )))
        if (is.na(type)) {
            stop(paste("Unknown method", method))
        } else if (type == 1) {
            interpfun <- interp.SPEC.LS
        }
        else if (type == 2) {
            interpfun <- interp.W.LK
        }
        multidiss.interp.SPEC(series, n, interpfun, integrate.divergenceW, alpha)
    } else {
        pairwise.diss( series, noindicesdiss(diss.SPEC.LLR), alpha, method, plot, n)
    }
}
multidiss.SPEC.GLK <- function(series, plot=FALSE) {
        multidiss.interp.SPEC(series, floor(length(series[[1]])/2), interp.SPEC.GLK, integrate.GLK)
}
multidiss.SPEC.ISD<- function(series, plot=FALSE,  n=length(series[[1]])) {
    if (n > 0) {
        multidiss.interp.SPEC(series, n, interp.SPEC.LOGLIKELIHOOD , integrate.ISD)
    } else {
        pairwise.diss( series, noindicesdiss(diss.SPEC.ISD), plot, n)
    }
}

    
diss.SPEC.GLK <- function(x,y, plot=FALSE ) {
    .ts.sanity.check(x, y)
    .check.equal.length.ts(x,y)
        
    interpx <- interp.SPEC.GLK(x, length(x)/2) #the value n is ignored noly used for compatibility with
    interpy <- interp.SPEC.GLK(x, length(x)/2) #multidiss.SPEC
    
    if (plot) {
        YksX <- log(spec.pgram(x,plot=FALSE)$spec)
        lambdasX <- spec.pgram(x,plot=FALSE)$freq
        YksY <- log(spec.pgram(y,plot=FALSE)$spec)
        lambdasY <- (spec.pgram(y,plot=FALSE)$freq)
        hX <- 0.93*dpill(lambdasX, YksX)
        hY <- 0.93*dpill(lambdasY, YksY)
        plotsmoothspec(lambdasX, YksX, lambdasY, YksY, .vectorized.lk.optim, hX, hY, length(x)/2)
    }
    #distance GLK
    integrate.GLK( interpx$x, interpx$y, interpy$y)
}



#distancia ISD
diss.SPEC.ISD <- function(x,y, plot=FALSE, n=length(x)) {
    .ts.sanity.check(x, y)
    .check.equal.length.ts(x,y)
 
    
    YksX <- log(spec.pgram(x,plot=FALSE)$spec)
    lambdasX <- spec.pgram(x,plot=FALSE)$freq
    YksY <- log(spec.pgram(y,plot=FALSE)$spec)
    lambdasY <- (spec.pgram(y,plot=FALSE)$freq)
    
    hX <- 0.93*dpill(lambdasX, YksX)
    hY <- 0.93*dpill(lambdasY, YksY)
    
    integraISDaux <- function(lambda) {
        (.vectorized.lk.optim(lambda, YksX, lambdasX  , hX) - .vectorized.lk.optim(lambda,YksY, lambdasY,  hY))**2
    }
    a <- 0
    if (n > 0) {
        a <- multidiss.interp.SPEC(list(x,y), n, interp.SPEC.LOGLIKELIHOOD , integrate.ISD)
    } else {
        tryCatch( {
            a <- integrate(integraISDaux, min(lambdasX), max(lambdasX) )$value
        }, error = function(e) {
            hX <- 2*hX
            hY <- 2*hY
            a <- integrate(integraISDaux, min(lambdasX), max(lambdasX) )$value
        })
    }
    if (plot) {
        plotsmoothspec(lambdasX, YksX, lambdasY, YksY, .vectorized.lk.optim, hX, hY, n)
    }  
    a
}


#########################################################
##############  distance CEPSTRAL  ######################
#########################################################


.calc.cepstral.coef <- function( ARx, h ) {
    CEPSTRALx <- 1:h
    CEPSTRALx[1] <- -ARx[1]
    if (length(ARx) >= 2) {
        for (i in 2:length(ARx)) {
            acum <- 0
            for (m in 1:(i-1) ) {
                acum <- acum + ( 1- m/i)*ARx[m]*CEPSTRALx[i-m]
            }
            CEPSTRALx[i] <- -ARx[i] + -acum
        }
    }
    
    if (h > length(ARx)) {
        for (i in (length(ARx)+1):h) {
            acum <- 0
            for (m in 1:length(ARx)) {
                acum <- acum + (1 - m/i)*ARx[m]*CEPSTRALx[i-m]
            }
            CEPSTRALx[i] <- -acum
        }
    }
    CEPSTRALx
    
}


cepstral <- function(x, h, order=NULL, seasonal, permissive) {
    ARx <- NULL
    
    SAR <- NULL
    if (is.null(order)) { #if order null, fit automatically
        ARx <- ar(x,order.max=min(length(x)-1,h))
    } else {
        if ((order[1]) < 1) stop("The arima order must have AR coefficients, they are used for the distance")
        arim <- arima(x, order, seasonal)
        ARx$ar <- arim$coef[1:order[1]]
        if (seasonal[[1]][1]>0) { #the seasonal part as in kalpakis
            SAR <- arim$coef[ (1 + order[1] + order[2]):(1 + order[1] + order[2]+ seasonal[[1]][1]) ]
        }
    }
    
    if (length(ARx$ar) < 1) {
        if (permissive) {
            warning("Cepstral distance, error on the selection of the AR order, 0 by AIC, forcing 1")
            ARx  <- ar(x, aic=FALSE, order.max=1)
        }
        if (length(ARx)<1) {
            stop("Could not find any AR coefficient")      
        }
    }
    
    cepst.SAR <- NULL #compute the sesonal part
    if (is.null(SAR)) {
        cepst.SAR <- rep.int(0, h)
    }
    else {
        cepst.SAR <- .calc.cepstral.coef(SAR, h)
    }
    
    rbind( .calc.cepstral.coef(ARx$ar, h), cepst.SAR )
    
}


diss.AR.LPC.CEPS <- function(x, y, k=50, order.x=NULL, order.y= NULL,
                    seasonal.x = list(order = c(0, 0, 0), period = NA),
                    seasonal.y = list(order = c(0, 0, 0), period = NA),
                    permissive=TRUE) {
    .ts.sanity.check(x, y)
    if (!is.list(seasonal.x) || !is.list(seasonal.y)) {
            stop("Invalid seasonal part")
    }
    cpx <- cepstral(x, k, order.x, seasonal.x, permissive)
    cpy <- cepstral(y, k, order.y, seasonal.y, permissive)
    as.numeric(dist(rbind(cpx[1,],cpy[1,])) + dist(rbind(cpx[2,],cpy[2,])))
}


#############################################################################
#################   Temporal Correlation Distance   #########################
#############################################################################

##CHOUAKRIA-DOUZAL

corrtemporder1 <- function (x, y) {
    p <- length(x)
    sum((x[2:p] - x[1:(p-1)]) * (y[2:p] - y[1:(p-1)])) / ( sqrt( sum((x[2:p] - x[1:(p-1)])^2) ) * sqrt( sum((y[2:p] - y[1:(p-1)])^2) ))
}

diss.CORT <- function( x, y, k=2, deltamethod="Euclid") {
    .ts.sanity.check(x, y)
    .check.equal.length.ts(x,y)
    corrt <- corrtemporder1(x,y)
    type <-  (pmatch(deltamethod, c("Euclid", "Frechet", "DTW")))
    typedist <- 0
    if (is.na(type)) {
        stop(paste("Unknown method", deltamethod))
    } else if (type == 1) {
        typedist <- as.numeric( dist(rbind(x,y)) )
    }
    else if (type == 2) {
        typedist <- diss.FRECHET(x,y)
    }
    else if (type == 3) {
        typedist <- dtw(x,y, dist.method="Manhattan", distance.only=T)$distance
    }
    
    (2/( 1+ exp(k*corrt)))*typedist
    
}



##################################################
######  maharaj clustering algorithm #############
##################################################
#input, distance matrix
pvalues.clust <- function(pvalues,significance) {
    distancias <- pvalues
    significacion <- significance
    distancias <- as.matrix(distancias)
    tam <- dim(distancias)[1]
    distancias
    combinaciones <- combn(tam,2)
    plandist <- 1:ncol(combinaciones)
    #create a vector with the distances
    for (i in 1:length(plandist)) {
        plandist[i] <- distancias[combinaciones[1,i], combinaciones[2,i]]
    }
    
    ord <- order(plandist, decreasing=TRUE)
    plandist <- plandist[ord]
    combinaciones <- combinaciones[,ord]
    
    is_in_setlist <- function( element, setlist) {
        for (i in 1:length(setlist)) {
            if (element %in% setlist[[i]]) {
                return(TRUE)
            }
        }
        return(FALSE)
    }
    
    find_which_set <- function ( element, setlist) {
        for (i in 1:length(setlist)) {
            if (element %in% setlist[[i]]) {
                return(i)
            }
        }
        return(0)
    }
    
    
    is_pvalue_of_element_less_than_significance_with_any_in_set <- function( element, setid, setlist, significance,distances) {
        sum(distances[element, setlist[[setid]]] < significance) > 0
    }
    
    add_to_set <- function ( element, setindex, setlist) {
        setlist[[setindex]] <- union(setlist[[setindex]], element)
        setlist
    }
    
    create_new_set <- function( element, setlist) {
        setlist[[length(setlist) + 1]] <- element
        setlist
    }
    
    is_all_series_already_in_a_cluster <- function(series, setlist) {
        already = TRUE
        for (i in series) {
            already <- already & (find_which_set(i,setlist) != 0)
        }
        already
    } 
    
    are_pvalues_of_all_pairs_across_clusters_greater_than_significance <- function( clusterone, clustertwo, setlist, significance, distances) {
        seriesone <- setlist[[clusterone]]
        seriestwo <- setlist[[clustertwo]]
        aregreater <- TRUE
        for (i in seriesone) {
            aregreater <- aregreater & (sum(( distances[i, seriestwo] > significance)) == length(seriestwo) )
        }
        for (i in seriestwo) {
            aregreater <- aregreater & (sum(( distances[i, seriesone] > significance)) == length(seriesone) )
        }
        aregreater
    }
    
    merge_sets <- function (setone, settwo, setlist) {
        setlist[[settwo]] <- c( setlist[[settwo]] ,setlist[[setone]])
        setlist[-setone]
    }
    
    if (plandist[1] < significacion) {
        grupos <- as.list(1:tam)
    } else {
        grupos <- list()
        grupos[[1]] <- combinaciones[,1]
        for (i in 2:length(combinaciones[1,])) {
            if (plandist[i] < significacion) { #p-value < significance = YES
                conj <- find_which_set( combinaciones[1,i], grupos )
                if (conj == 0) {
                    grupos <- create_new_set(combinaciones[1,i], grupos)
                }
                conj <- find_which_set( combinaciones[2,i], grupos )
                if (conj == 0) {
                    grupos <- create_new_set(combinaciones[2,i], grupos)
                } 
                #each remaining serie to its own cluster
                for (j in i:length(combinaciones[1,])) {
                    conj <- find_which_set( combinaciones[1,j], grupos )
                    if (conj == 0) {
                        grupos <- create_new_set(combinaciones[1,j], grupos)
                    }
                    conj <- find_which_set( combinaciones[2,j], grupos )
                    if (conj == 0) {
                        grupos <- create_new_set(combinaciones[2,j], grupos)
                    } 
                }
                break;
            } else { #p-value < significance = NO
                conj <- is_in_setlist( combinaciones[1,i], grupos )
                conj <- conj + is_in_setlist( combinaciones[2,i], grupos )
                if (conj < 2) {  #is each(ALL) series already in a cluster = NO
                    conj <- find_which_set( combinaciones[1,i], grupos )
                    if (conj > 0) { #one of the series in the pair already in a cluster = YES (x)
                        if (is_pvalue_of_element_less_than_significance_with_any_in_set( combinaciones[2,i], conj, grupos, significacion, distancias  )    ) {
                            grupos <- create_new_set(combinaciones[2,i], grupos)
                        } else {
                            grupos <- add_to_set( combinaciones[2,i], conj, grupos )
                        }
                    } else {
                        conj <- find_which_set( combinaciones[2,i], grupos )
                        if (conj > 0) {#one of the series in the pair already in a cluster = YES (y)
                            if (is_pvalue_of_element_less_than_significance_with_any_in_set( combinaciones[1,i], conj, grupos, significacion, distancias  )    ) {
                                grupos <- create_new_set(combinaciones[1,i], grupos)
                            } else {
                               grupos <- add_to_set( combinaciones[1,i], conj, grupos )
                            }
                        } else {#one of the series in the pair already in a cluster = NO
                            grupos <- create_new_set( combinaciones[1,i], grupos ) # create a new cluster with the two series
                            setid <- find_which_set( combinaciones[1,i], grupos )
                            grupos <- add_to_set( combinaciones[2,i], setid, grupos )
                        }
                    } 
                } else {  #is each series already in a cluster = YES
                    conj1 <- find_which_set( combinaciones[1,i], grupos )
                    conj2 <- find_which_set( combinaciones[2,i], grupos )
                    if (are_pvalues_of_all_pairs_across_clusters_greater_than_significance(conj1, conj2, grupos, significacion, distancias)) {
                        #merge
                        grupos <- merge_sets( conj1, conj2, grupos)
                    } else {
                        #do nothing
                    }
                }
            }
        }  
    }
    
    result <- 1:tam
    for ( i in 1:length(grupos) ) {
        for (j in grupos[[i]]) {
            result[j] <- i
        }
    }
    result
}


#check if the maharaj algorithm is equivalent to hierarchical clustering with 
# a cut level
testIgualdadMaharajHCLUST <- function( distancias, pvalor) {
    clusterp <- pvalues.clust(as.dist(distancias), pvalor)
    clusterh <- (hclust(as.dist(-distancias), method="single"))
    print(clusterp)
    plot(clusterh)
}


####################################################################
############# Feature Extraction Based on Wavelet Transform ########
####################################################################

wavelet.feature.extraction <- function(series) {
    
    calcEnergies <- function(wavdecomp) {
        level <- length(wavdecomp$data) -1
        energyD <- rep.int(0,level)
        energyD[1] <- sum(wavdecomp$data[[1]]**2)
        for ( i in 1:level) {
            energyD[i] <- sum((wavdecomp$data[[i]])**2)
        }
        return (energyD)
    }
    
    #fill with zeroes
    max.level <- as.integer(ceiling(logb(length(series[1,]),base=2)))
    true.level <- as.integer(floor(logb(length(series[1,]),base=2)))
    if (max.level != true.level) {
        npad <- 2**max.level - length(series[1,])
        for (i in 1:npad) {
            series <- cbind(series, rep(0, nrow(series)))
        }
    }
    
    energies <- matrix(0, nrow=nrow(series), ncol = max.level)
    for ( i in 1:nrow(series) ) {
        wavdecomp <- wmtsa::wavDWT(series[i,], n.levels=max.level, wavelet="haar")
        energies[i,] <- calcEnergies(wavdecomp)
    }
    
    sumEnergies <- colSums(energies)
    
    final_level <- max.level
    for (i in 1:(max.level-1)) {
        if (sumEnergies[i] < sumEnergies[i+1]) {
            final_level <- i
            break
        }
    }
    wavdecomp <- wavDWT(series[1,], n.levels=final_level, wavelet="haar")
    out.series <- wavdecomp$data[[final_level+1]]
    for (i in 2:nrow(series)) {
        wavdecomp <- wavDWT(series[i,], n.levels=final_level, wavelet="haar")
        out.series <- rbind(out.series, wavdecomp$data[[final_level+1]])
    }
    out.series
}

diss.DWT <- function(series) {
    if ( length(dim(series)) == 2 ) {
        if ( dim(series)[1] < 2 ) {
            stop( "diss.DWT needs a minimum of 2 series to compute the distance, incorrect amount provided" )
        }
    } else {
        stop( "diss.DWT series matrix with incorrect dimensions" ) 
    }
    wt <- wavelet.feature.extraction( series )
    dist(wt)
}

##########################################################
############  CORRELATION BASED DISTANCES ################
##########################################################



diss.COR <- function(x, y, beta = NULL) {
    .ts.sanity.check(x, y)
    correl <- cor(x,y)
    if (is.null(beta)) {
        sqrt(2*(1- correl))
    } else {
        if (beta<0) {
            stop("beta must be greater than 0")
        }
        sqrt( ((1-correl)/(1+correl ))**beta )    
    }
}


####################################################################################
########################## CDM KEOGH 2004 ##########################################
################## Compression based data mining of sequential data ################
####################################################################################

#common part of compression methods,
#calculate the sizes of the compressed series and of their concatenation
.compression.lengths <- function(x, y, type) {       
    methods <- type
    type = match.arg(type, c("gzip", "bzip2", "xz", "min"))
    if (type == "min") { #choose the best compression method of the three 
        methods <- c("gzip", "bzip2", "xz")
    }
    xy <- as.character(c(x,y))
    x <- as.character(x)
    y <- as.character(y)
    cxym <- sapply( methods, function(m) { length( memCompress(xy, type=m) )})
    cxy <- min(cxym)
    cx <- min(sapply( methods, function(m) { length( memCompress(x, type=m) )}))
    cy <- min(sapply( methods, function(m) { length( memCompress(y, type=m) )}))
    list(cx=cx, cy=cy, cxy=cxy)
}

#length of the compressed concatenated series / sum lengths of the compressed series
diss.CDM <- function(x, y, type="min") {
    .ts.sanity.check(x, y)
    comp <- .compression.lengths(x,y, type)    
    comp$cxy / (comp$cx + comp$cy)
}


###################################################################################
####### Clustering by Compression (2005), Cilibrasi, R., Vitanyi, P.M.B.,  ########
######## Normalized Compression Distance ##########################################
###################################################################################
diss.NCD <- function(x,y, type="min") {
    .ts.sanity.check(x,y)
    comp <- .compression.lengths(x,y, type)  
    (comp$cxy - min(comp$cx,comp$cy)) / max(comp$cx, comp$cy)
}

##############################################################
################# COMPLEXITY BASED DISTANCE ##################
#BatistaWangKeogh_2011_SDM_A complexity invariant distance measure for TS#
##############################################################

diss.CID = function(x, y) {
    .ts.sanity.check(x, y)
    .check.equal.length.ts(x,y)
    CED.x = sqrt( sum( diff(x)^2) ) #complexities of the series
    CED.y = sqrt( sum( diff(y)^2) )
    CF = max(CED.x, CED.y) / min(CED.x, CED.y) #complexity correction factor
    CF * dist(rbind(x,y))
}



#######################################################################
############### CLUSTER SIMILARITY INDEX ##############################
#######################################################################

#gravrilov similarity ratio, used in kalpakis
Sim <- function(Gi, Sj, G, S) {
    2* (sum ( (G==Gi) & (S==Sj) ) ) / ( sum(G==Gi) + sum(S==Sj))
}

cluster.evaluation <- function(G,S) {
    if (length(G) != length(S)) {
        stop("Different amount of elements between cluster solutions")
    }
    if (any(is.na(G)) || any(is.na(S))) {
        stop("NA in the cluster solutions")
    }
    acum <- 0
    gclust <- unique(G)
    sclust <- unique(S)
    for (i in gclust) {
        maxS <- 0
        for (j in sclust) {
            ms <- Sim(i, j, G, S)
            if (ms > maxS) {
                maxS <- ms
            }
        }
        acum <- acum + maxS
    }
    acum/(length(gclust))
}


###################################################################
########################### DISS WRAPPER ##########################
###################################################################

#series is a list with time series
#dissfun is a function that takes the list of series and their indices, and extra paramters
#... parameters for dissfun
pairwise.diss <- function(series, dissfun, ...) {    
    n <- length(series)
    distances <- matrix(0, n, n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            tryCatch( {
                d <- dissfun( series, i, j, ...)
                distances[i,j] <- d
                distances[j,i] <- d
            }, error = function (e) {
                stop( paste("Applying diss, series (",i,",",j,") produced the following error: ", e) )
            })
        }
    }
    as.dist((distances))
}

noindicesdiss <- function( fun ) {
    function(series, i, j, ...) {
        fun(series[[i]], series[[j]], ...)
    }
}

diss <- function(SERIES, METHOD, ...) {
    if (!is.matrix(SERIES) && !is.list(SERIES) && !is.mts(SERIES)) {
        stop("list, mts, matrix or data.frame object is required for SERIES ")
    }
    mat.ser <- SERIES
    if (is.mts(SERIES)) {
        SERIES <- t( as.matrix(SERIES))
    }
    
    if (!is.list(SERIES)) {
        tmpser <- SERIES
        SERIES <- list()
        for (i in 1:nrow(tmpser)) {
            SERIES[[i]] <- tmpser[i,]
        }
        names(SERIES) <- rownames(tmpser)
    }
    
    #common check for input parameters
    if (any(is.na(unlist(SERIES)))) {
        stop("NA in the series")
    }
    if (length(SERIES) < 2) {
        stop("Only one series provided")
    }
    
    list.to.matrix <- function(series) {
        n <- length(series)
        k <- length(series[[n]])
        mat.ser <- matrix(0, n, k)
        for (i in 1:n) {
            if ( length( series[[i]]) != k ) {
                stop("diss method requires same length series")
            }
            mat.ser[i,] <- series[[i]]
        }
        rownames(mat.ser) <- names(series)
        mat.ser
    }

    out.dist <- NULL
    
    METHODS = c("ACF", "PACF", "AR.MAH", "AR.PIC", "AR.LPC.CEPS", "PER", "INT.PER", "COR", "CORT", "DWT",
                "PDC", "PRED", "MINDIST.SAX", "SPEC.LLR", "SPEC.GLK", "SPEC.ISD", "CDM", "CID", "NCD", "DTWARP", "FRECHET",
                "EUCL")
    diss.method = match.arg(METHOD, METHODS)
    #get the statistic of the MAHARAJ dissimilarity
    diss.AR.MAH.STAT <- function(x,y, ...) {
        diss.AR.MAH(x,y,...)$statistic
    }
    diss.AR.MAH.PVAL <- function(x,y, ...) {
        diss.AR.MAH(x,y,...)$p_value
    }
    diss.fun <- switch(diss.method,
        ACF = diss.ACF,
        PACF = diss.PACF,
        AR.MAH = diss.AR.MAH,
        AR.PIC = diss.AR.PIC,
        AR.LPC.CEPS = diss.AR.LPC.CEPS,
        PER = diss.PER,
        INT.PER = diss.INT.PER,
        COR = diss.COR,
        CORT = diss.CORT,
        DWT = diss.DWT,
        PRED = multidiss.PRED,
        SPEC.LLR = diss.SPEC.LLR,
        SPEC.GLK = diss.SPEC.GLK,
        SPEC.SD = diss.SPEC.ISD,
        MINDIST.SAX = diss.MINDIST.SAX,
        PDC = pdcDist,
        CDM = diss.CDM,
        CID = diss.CID,
        NCD = diss.NCD,
        DTWARP = diss.DTWARP,
        FRECHET = diss.FRECHET,
        EUCL = diss.EUCL)
    
    if (diss.method == "DWT") { #diss dwt is not a pairwise diss, we cannot use proxy::dist
        out.dist <- diss.DWT( list.to.matrix(SERIES) )
    } else if (diss.method == "AR.PIC") {
        multi.PIC <- function(SERIES, i, j, order=NULL, permissive=TRUE, order.x=NULL, order.y=NULL) {
            if (!is.null(order.x) || !is.null(order.y) ) {
                stop("AR.PIC from the diss wrapper function must be called using 'order' argument, not with
                     order.x and order.y, see diss.AR.PIC help page")
            }
             diss.AR.PIC(SERIES[[i]], SERIES[[j]], order[i,], order[j,], permissive)
        }
        out.dist <- pairwise.diss(SERIES, multi.PIC, ...)
    } else if (diss.method == "AR.LPC.CEPS") {
        multi.CEPS <- function(series, i, j, k=50, order=NULL, seasonal=NULL, permissive=TRUE,
                        order.x=NULL, order.y=NULL, seasonal.x=NULL, seasonal.y=NULL) { #arguments to inform incorrect usage
            if (!is.null(order.x) || !is.null(order.y) || !is.null(seasonal.x) || !is.null(seasonal.y)) {
                stop("AR.LPC.CEPS from the diss wrapper function must be called using 'order' and 'seasonal' arguments, not with
                     order.x, order.y, seasonal.x or seasonal.y arguments, see diss.AR.LPC.CEPS help page")
            }
            if (is.null(seasonal)) {
                seasonal[[i]] <-  list(order=c(0,0,0), period=NA)
                seasonal[[j]] <-  list(order=c(0,0,0), period=NA)
            }
            distance <- diss.AR.LPC.CEPS(series[[i]], series[[j]], k, order[i,], order[j,], seasonal[[i]], seasonal[[j]] )
        }
        out.dist <- pairwise.diss(SERIES, multi.CEPS, ...)
    } else if (diss.method == "PRED") {
        return( multidiss.PRED(SERIES, ...) ) #TODO proper names
    } else if (diss.method == "AR.MAH") {

        statistic = pairwise.diss( SERIES, noindicesdiss(diss.AR.MAH.STAT), ...)
        p_value = pairwise.diss( SERIES, noindicesdiss(diss.AR.MAH.PVAL), ...)
        return( list(statistic=statistic, p_value=p_value) ) #TODO proper naming of the output
    } else if (diss.method == "PDC") {
        out.dist <- pdcDist( t( list.to.matrix(SERIES) ), ...)
    } else if (diss.method == "SPEC.LLR") { #for performance reasons, we must call these in a different way
        out.dist <- multidiss.SPEC.LLR(SERIES, ...)
    } else if (diss.method == "SPEC.GLK") { #for performance reasons, we must call these in a different way
        out.dist <- multidiss.SPEC.GLK(SERIES, ...)
    } else if (diss.method == "SPEC.ISD") { #for performance reasons, we must call these in a different way
        out.dist <- multidiss.SPEC.ISD(SERIES, ...)
    } else {
        out.dist <- pairwise.diss( SERIES, noindicesdiss(diss.fun), ...)
    }
    out.dist <- as.dist(out.dist)
    names(out.dist) <- names(SERIES)
    out.dist
}

loo1nn.cv <- function(d, G) {
    d <- as.matrix(d)
    diag(d) <- max(d) + 1 #distance with self is not included
    nearest.set <- apply ( d, 1, function(x) { #handle ties
                                    which( x == min(x) )
                                }
                         )
    
    #mode auxiliary function
    #based on http://stackoverflow.com/questions/2547402/standard-library-function-in-r-for-finding-the-mode
    #with random tie breaking
    my.mode <- function(x) {
        ux <- unique(G[x])
        if (length(ux) > 1) {
            ux <- sample(ux)
        }
        prop <- tabulate(match(G[x], ux))
        if (sum(prop == max(prop))>1) {
            warning("There were ties on the voting, selecting one at random")
        }
        ux[which.max(prop)]
    }
    
    #apply the mode, depending on whether apply obtained different amount of ties (list), the same amount(matrix) or no ties(vector) 
    nearest <- switch( class(nearest.set), matrix = apply(nearest.set, 2, my.mode),
                       numeric = apply(as.matrix(nearest.set), 1, my.mode),
                       integer = apply(as.matrix(nearest.set), 1, my.mode),
                       list = sapply(nearest.set, my.mode)
    )
    sum(nearest == G)/length(G)
}

#wrappers for easier discovery of the available functions
diss.EUCL <- function(x, y) {
    dist(rbind(x,y))
}
diss.DTWARP <- function(x,y,...) {
    dtw(x,y, ...)$distance
}
diss.FRECHET <- function(x,y,...) {
    abscissex = 1:length(x)
    abscissey = 1:length(y)
    distFrechet(abscissex,x,abscissey,y, ...)
}
diss.PDC <- function(x,y, ...) {
    pdcDist(cbind(x,y), ...)
}
############################################################################
#######################   OLD STUFF (UNUSED)   #############################
############################################################################



#distancia basada en correlaciones cruzadas
distanciaCorCruLagK = function(x,y,k) {
    muX = mean(x)
    muY = mean(y)
    TT = length(x)
    sum( (x[1:(TT-k)] - muX)*( y[(1+k):TT] - muY)) / (sqrt( sum((x - muX)**2))*sqrt( sum((y - muY)**2)))
}

#distanciaCorCruLagK(x,y,18)

distanciaCorCruTotal = function(x,y) {
    k = length(x)-1
    denom = 0
    for (i in 1:k) {
        denom = denom + distanciaCorCruLagK(x,y,i)
        print(denom)
    }
    
    sqrt( (1 - distanciaCorCruLagK(x,y,0))/denom  )
    denom
}

#distanciaCorCruTotal(x,y)
