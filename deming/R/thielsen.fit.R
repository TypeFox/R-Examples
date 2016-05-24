theilsen.fit <- function(x, y, wt, symmetric, conf, eps) {
    # Form all pairs of observations
    pdiff <- function(x, fun='-') {
        # give all paired distances
        n <- length(x)
        indx1 <- rep(1:(n-1), (n-1):1)
        indx2 <- unlist(lapply(2:n, function(x) x:n))
        as.vector((get(fun))(x[indx2], x[indx1]))  #no labels
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
    if (all(wt==1)) weighted <- FALSE
    else {
        ww <- pdiff(wt, "*")
        weighted <- TRUE
    }

    # Remove ties
    xtie <- (abs(xx) < eps)
    if (symmetric) tied <- (xtie & abs(yy)<eps)
    else    tied <- xtie
    if (any (tied)) {
        xx <- xx[!tied]
        yy <- yy[!tied]
        if (weighted) ww <- ww[!tied]
        xtie <- xtie[!tied]   #symmetric option will note tied x values
    }

    theta <- atan(yy/xx)  # angle
    if (!weighted) {  #the most common case
        if (symmetric) {
            # Find that ray with angle 0 <= gamma <= pi/2 such that 1/2 of the
            #  thetas lie in the right angle formed by gamma and gamma-pi/2
            # If the overall slope is negative first reflect the points about
            #  the x=0 axis.
            slope <- tan(median(theta[!xtie]))
            if (slope <0) theta <- -theta
            theta[xtie] <- - pi/2  # a special case
            npair <- length(theta)
            # sweep the angle through the cloud, keeping track of how
            #  many points are within it.  There may be multiple solutions
            #  so approx() is not appropriate
            # angle2 = angle at which each point changes allegiance
            angle2 <- ifelse(theta >=0, theta, theta + pi/2)
            indx <- order(angle2)
            ss <- (ifelse(theta <0, -1, 1))[indx]
            angle2 <- angle2[indx]
            init <- sum(theta < 0)
            if (min(angle2) > 0) {
                angle2 <- c(0, angle2)
                ss <- c(0, ss)
            }
            # inside = number inside the right angle, with those on the line
            #  counting as 1/2, minus npair/2
            inside <- as.vector((init + cumsum(ss) - ss/2)) - npair/2  #no names
            cross <- which(abs(diff(sign(inside))) ==2) #zero crossings
            z <- cross+1
            zeros <- (angle2[cross]*inside[z] - angle2[z]*inside[cross]) /
                      (inside[z] - inside[cross])
            zeros <- c(angle2[inside==0], zeros)
            if (length(zeros) ==1) slope <- tan(zeros) * ifelse(slope<0, -1,1)
            else {
                # Choose the solution with minimal rotated residual
                newy <- cbind(y,x) %*% rbind(cos(zeros), -sin(zeros))
                asum <- apply(newy, 2, mad)
                slope <- tan(zeros[order(asum)][1])
            }
            return(list(coefficients=c(median(y - x*slope), slope),
                        angle = zeros))
        }
        else if (conf > 0) {
            slope <- tan(median(theta))
            # Compute the standard Theil-Sen confidence interval
            # Never done for the symmmetric=TRUE case
            n1 <- length(x)  #the n of the data sample
            npair <- length(theta)
            v = sqrt(n1* (n1-1) *(2*n1 +5)/ 18) # Sen, equation 2.6
            tiecount <- as.vector(table(x))
            v <- v - sum(tiecount*(tiecount-1)* (2*tiecount +5))/18
            z <- -qnorm((1-conf)/2)
            dist <- ceiling(v * z/2) # this many points above and below
            ci <- approx(1:npair -.5, sort(theta), npair/2 +c(-dist, dist))$y
            ci <- tan(ci)
            ci <- matrix(c(median(y - x*ci[1]), median(y- x*ci[2]), ci),
                         byrow=TRUE, ncol=2)
          return(list(coefficients=c(median(y - x*slope), slope), ci=ci))
        }
        else{
            slope <- tan(median(theta))
            return(list(coefficients=c(median(y - x*slope), slope)))
        }
    }

    # The weighted case, no confidence intervals
    ww <- ww/mean(ww)  #weights now sum to npair
    slope <- tan(wtmed(theta, ww))
    if (symmetric) {
        if (slope <0) theta <- - theta
        theta[xtie] <- - pi/2  # a special case
        npair <- length(theta)
        angle2 <- ifelse(theta >=0, theta, theta + pi/2)
        indx <- order(angle2)
        ss <- (ifelse(theta<0, -ww, ww))[indx]
        angle2 <- angle2[indx]
        if (min(angle2)>0) {
            angle2 <- c(0, angle2)
            ss <- c(0,ss)
        }
        init <- sum(ww[theta <0])
        inside <- as.vector((init + cumsum(ss) - ss/2) - npair/2)
        cross <- which(abs(diff(sign(inside))) ==2) #zero crossings     
        z <- cross+1
        zeros <- (angle2[cross]*inside[z] - angle2[z]*inside[cross]) /
                      (inside[z] - inside[cross])
        zeros <- c(angle2[inside==0], zeros)
        if (length(zeros)==1) slope <- tan(zeros * ifelse(slope<0, -1, 1))
        else {
            # Choose the solution with minimal rotated residual
            newy <- cbind(y,x) %*% rbind(cos(zeros), -sin(zeros))
            asum <- apply(newy, 2, function(z)
                          wtmed(abs(z-wtmed(z, wt)),wt))
            slope <- tan(zeros[order(asum)][1])
        }
        list(coefficients=c(wtmed(y - x*slope, wt), slope),
             angle = zeros)
    }
    else list(coefficients=c(wtmed(y- x*slope, wt), slope))
}
