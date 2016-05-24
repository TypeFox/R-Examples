pbreg.fit<- function(x, y, wt, method, conf, eps) {
    pdiff <- function(x, fun='-') {
        # give all paired distances
        n <- length(x)
        indx1 <- rep(1:(n-1), (n-1):1)
        indx2 <- unlist(lapply(2:n, function(x) x:n))
        as.vector((get(fun))(x[indx2], x[indx1]))  #drop any names
    }

    # weighted median.  Essentially a CDF of the data, then connect
    #  midpoints of the rise
    wtmed <- function(y, wt) {
        indx <- order(y)
        ww <- wt[indx]
        approx(cumsum(ww) - ww/2, y[indx], sum(ww)/2)$y
    }

    xx <- pdiff(x)
    yy <- pdiff(y)
    if (all(wt==wt[1])) weighted <- FALSE
    else {
        weighted <- TRUE
        ww <- pdiff(wt, "*")
    }
    
    uninformative <- (abs(xx) < eps & abs(yy) < eps)
    if (any(uninformative)) {
        xx <- xx[!uninformative]
        yy <- yy[!uninformative]
        if (weighted) ww <- ww[!uninformative]
    }

    # change to polar coordinates
    # negative values are clockwise from the horizontal axis
    theta <- atan(yy/xx)  #ranges from -pi/2 to pi/2
    if (method==1) {
        #Passing-Bablock, symmetric around 45 degree line
        theta <- ifelse(theta < -pi/4, theta + pi, theta)
        keep <- (abs(xx + yy) > eps)  #not on the -45 degree (-pi/4) line
    }
    else if (method==2) {
        #Passing-Bablock method 2
        below <- (theta <0 & abs(xx) > eps) #don't count tied x values
        if (any(below)) {
            if (weighted) temp <- wtmed(theta[below], ww[below])
            else  temp <- median(theta[below])
        }
        else temp <- -1  #dummy value for the rare case of monotone data
        theta <- ifelse(theta < temp, theta + pi, theta)
        keep <- (abs(xx*cos(temp) + yy*sin(temp)) > eps)
    } else {
        # The Passing-Bablock "scissors" estimator
        theta = abs(theta)
        keep <- rep(TRUE, length(theta))
    }
    theta <- theta[keep]
    if (!weighted) {  # the usual case
        slope <- tan(median(theta))
        if (conf>0) {
            # Compute the standard Theil-Sen confidence interval
            n1 <- length(x)  #the n of the data sample
            npair <- length(theta)
            v = sqrt(n1* (n1-1) *(2*n1 +5)/ 18) # Sen, equation 2.6
            tiecount <- as.vector(table(x))
            v <- v - sum(tiecount*(tiecount-1)* (2*tiecount +5))/18
            z <- -qnorm((1-conf)/2)
            dist <- ceiling(v * z/2) # this many points above and below
            ci <- approx(1:npair -.5, sort(theta), npair/2 +c(-dist, dist))$y
            ci <- tan(ci)
            ci <- matrix(c(median(y - x*ci[2]), median(y- x*ci[1]), ci),
                         byrow=TRUE, ncol=2)
          list(coefficients=c(median(y - x*slope), slope), ci=ci)
        }
        else list(coefficients=c(median(y - x*slope), slope))
    }
    else {                  
        # The less common weighted case
        ww <- ww[keep]
        slope <- wtmed(theta, ww)
        list(coefficients=c(wtmed(y - x*slope, wt), slope))
    }     
}
