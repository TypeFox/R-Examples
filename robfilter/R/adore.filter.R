adore.filter <- function(y, 
                         p.test = 15, 
                         minNonNAs = 5,
                         min.width = 10, 
                         max.width = 200, 
                         width.search="geometric",
                         rtr=2, 
                         extrapolate=FALSE, 
                         calc.qn = FALSE, 
                         sign.level=0.1
                        ) {

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Validity check
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    if(missing(y)){
        stop("The input data vector is missing with no default.\n")
    } else {
        tseries <- y   
        rm(y)
    }
    
    if(!is.numeric(tseries))
        stop("The data vector must be numeric.\n")
        
    if(!is.numeric(min.width))
        stop("The minimal window width (min.width) must be numeric.\n")

    if(min.width < 5)
        stop("The minimal window width allowed is min.width=5.\n")
    
    if(!(rtr %in% c(0,1,2)))
        stop("The parameter rtr must be either 0, 1 or 2.\n
              0: no restriction of the level estimate to the observational range\n
              1: restriction of the level estimate to the observational range within the current window\n
              2: restriction of the level estimate to the observational range of the most recent observation within the current window")

    N <- length(tseries) # length of the time series
    if(N < min.width)
        stop("The data vector must be longer than min.width.\n")

    if(!is.numeric(minNonNAs))
        stop("minNonNAs must be numeric.\n") 
                
    if(!is.numeric(max.width) && !(is.logical(max.width) && !max.width))
        stop("The maximal window width (max.width) must be either numeric or FALSE.\n")

    if (!(is.logical(max.width) && !max.width)){ 
      if(min.width > max.width)
        stop("The maximal window width (max.width)must be >= min.width or FALSE.\n")
    }

    if(! (p.test %in% c(0.25, 0.3, 0.5, 5:120)))
        stop("p.test must be either 0.25, 0.3 or 0.5 as a fraction of
             observations for each window or a fixed integer in 5,...,120.\n")

    if (p.test < 1){
      if ((minNonNAs < 5) || (minNonNAs > floor(min.width/2)))
        stop("minNonNAs must be numeric and 5 <= minNonNAs <= floor(min.width/2).\n") 
      nI.start <- floor(min.width/2)
    } else {
      if ((minNonNAs < 5) || (minNonNAs > p.test))
        stop("minNonNAs must be numeric and 5 <= minNonNAs <= p.test.\n") 
        nI.start <- min(p.test, floor(min.width/2))
    }

    if( !(width.search %in% c("linear","binary","geometric","set.to.min"))){
        stop("width.search must be one of 'linear', 'binary', 'geometric' or 'set.to.min'.\n")
    }
    
    if(sign.level <= 0 | sign.level > 0.5)
        stop("'sign.level' must be a value in (0,0.5)")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Starting window width
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    start <- min.width    
    n.non.missing <- ifelse( (nI.start < minNonNAs), nI.start, minNonNAs )
    # Enlarging the first window width if there are too many missing values in the first window:
    if(sum(! is.na(tseries[(min.width-nI.start+1):min.width])) < n.non.missing) {
        warning("The first time window must contain at least ", n.non.missing, 
                " observations at the most recent time points.\n")

        n.s <- c(rep(NA,min.width-1),min.width:N)
        if(p.test >= 5){
            nI.s <- apply(cbind(rep(p.test,N),floor(n.s/2)),1,min)
        } else {
            nI.s <- apply(cbind(floor(p.test*n.s), rep(floor(min.width/2),N)),1,max)
        }
        
        while( sum(!is.na(tseries[(start-nI.s[start]+1):start])) < ifelse( (nI.s[start] < minNonNAs), nI.s[start], minNonNAs ) ){
            start <- start+1
            if(start > N) stop("The time series must contain >=", minNonNAs, " non-NAs.\n")
        }

        warning("The initial window width has been enlarged to ", start, ".\n")
    }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Internal Functions
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Loading of critical values for the test of adequacy of the current fit in each time window
    #data(critvals)
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Function to get the suitable critical value at each time
    
    get.critval <- function(n, nI) {
        if( sign.level==0.1 & n <= 600  & nI <= 61 & n >=11 & nI >= 5){
            return(critvals[n, nI])
        } else {
            k <- n %% 2
            return(2 * qhyper(1 - sign.level/2, (n - k)/2, (n - k)/2, nI) - nI)
        }
    }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Function to calculate the test statistic
    get.T <- function(res) {
        return(abs(sum(sign(res), na.rm = TRUE)))
    }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Function to calculate the value of nI within each window (must be >=5)
    get.nI <- function(p.test, n) {
        if(p.test < 1)
            return(floor(p.test * n))
        if(p.test > n/2)
                return(floor(n / 2))
        return(p.test)
    }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Function to calculate the residuals, also considering rounding errors
    # -> tolerance (tol) depending on the System, on a 64 bit: approx. 2.220446e-16^0.5
    get.res <- function(y, mu, beta, n, tol = .Machine$double.eps^0.5) {
        res <- y - (mu + beta * (-n + 1):0)
        return(replace(res, abs(res) < tol, 0))
    }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Correction in case of mu[t] > max(tseries[I.Win/I.t]) || mu[t] < min(tseries[I.Win/I.t])
    # Restricts possible signal estimations to observational range
    # Note: The last value in the vector mu.list[[t]] still corresponds to the signal estimation
    #       without the 'restrict to observational range' rule and hence might differ from mu[t]
    restrict.mu <- function(mut, y, I.t, rtr){
        if(rtr == 1){
            mu.range <- range(y, na.rm = TRUE)
        }
        if(rtr == 2){
            mu.range <- range(y[I.t], na.rm=TRUE)
        }  
        if(mut < mu.range[1]){
           mut <- mu.range[1]
        } 
        if(mut > mu.range[2]){ 
           mut <- mu.range[2] 
        }
        return(mut)  
    }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Internal function for evaluation of the adequate window width + corresponding RM regression fit in this window
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
get.fit.and.n <- function(y=y, M=M){
    # Internal variables:
    n      <- length(y)
    H      <- NULL # vector with the test decision of each iteration step
    i      <- 0    # initialising iteration step counter
    n.list <- n    # vector of window widths in each iteration step
    n.low  <- min.width # window width 'boundaries'
    n.up   <- n
    y0     <- y    # initial observations in window
    M0     <- M    # initial matrix of all pairwise observational slopes in the window

    all.levels.t <- NULL
    all.slopes.t <- NULL
    all.widths.t <- NULL

    repeat { # loop for window width adaptation
        i   <- i+1 # iteration step counter
        n.list[i] <- n
        nI  <- get.nI(p.test, n) # number of residuals for calculating the test statistic
        I.t <- (n - nI + 1):n

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Repeated Median estimation for more than minNonNAs observations are present at recent nI times
# else: NAs are returned for level and slope and the window width is not enlarged
        if( sum(! is.na(y[I.t])) < ifelse( (nI < minNonNAs), nI, minNonNAs ) ){
            level.t <- NA
            slope.t <- NA
            nt      <- n
            all.levels.t <- NA
            all.slopes.t <- NA
            all.widths.t <- NA
            break
        }

        betas <- apply(M, 1, median, na.rm = TRUE)
        if(any(! is.na(betas))) {
            slope.t <- median(betas, na.rm = TRUE)
            level.t <- median(y - slope.t * (-n + 1):0, na.rm = TRUE)

            all.levels.t[i] <- level.t
            all.slopes.t[i] <- slope.t
            all.widths.t[i] <- n
        }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Calculation of residuals
        res <- get.res(y, level.t, slope.t, n)
       
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Test decision
         reject.H0 <- (get.T(res[I.t]) > get.critval(n, nI))
         H[i]      <- reject.H0 # Saving the test decision (0: do not reject H0, 1: reject H0)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Choice of new window width with geometric search
        if(width.search=="geometric"){

            if( reject.H0 ){ # H0 rejected:

                if(n == min.width){
                    nt <- n # final window width = min.width
                    break   # H0: "Fit is good" is rejected but window width cannot be reduced anymore
                }
                # 'new' window width for the next iteration (binary search)
                n.up <- n
                if( any(!H) ){
                    n <- ceiling( mean(c(n.low,n.up)) ) 
                } else {
                    n    <- max(n - 2^(i-1), min.width)
                }

            } else { # H0 not rejected:

                if(i==1){
                    nt <- n # current width= final window width
                    break   # because H0: "Fit is good" cannot be rejected
                }
                # 'new' window width for the next iteration (binary search)
                n.low <- n
                n <- ceiling( mean(c(n.low,n.up)) ) 

            } # end of H0 not rejected

            # Final window width found:
            if( (n%in% n.list) & ((n==n.low) | (n==n.up)) ){

                max.n   <- max(which(!H))
                nt      <- all.widths.t[max.n]  # final window width
                n       <- nt
                slope.t <- all.slopes.t[max.n]  # final slope estimate
                level.t <- all.levels.t[max.n]  # final signal estimate
                w <- c((length(y0)-nt+1):length(y0))
                y <- y0[w]    # data vector for next step
                M <- M0[w,w]  # matrix with slopes in next step
                break
            } 
            
            # internals for next iteration
            w <- c((length(y0)-n+1):length(y0))
            y <- y0[w]   
            M <- M0[w,w] 

        } # end of geometric search

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Choice of new window width with linear one-step search
        if(width.search=="linear"){
            if( (n == min.width) || !reject.H0 ){
                nt   <- n
                break  # H0: "Fit is good" cannot be rejected or window width cannot be reduced
            } else {   # H0: "Fit is good" is rejected
                n     <- n - 1     # -> reduce window width by 1 and re-estimate level and slope in smaller window
                y     <- y[-1]     # data vector in next iteration
                M     <- M[-1, -1] # matrix with slopes in next iteration
            }
        } # end of linear one-step search

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Choice of new window width with binary search
        if(width.search=="binary"){
             if(i ==1){
                if( (n != min.width) & reject.H0 ){
                     n <- min.width  # fit next time in window with minimal window width
                } else  {
                    nt <- n
                    break        # fit with current window width is appropriate
                }
            # further iterations:
            } else {
                if(reject.H0){
                    if(n == min.width){
                        nt   <- n
                        break  # Although H0: "Fit is good" is rejected, the window width cannot be reduced anymore
                    }
                    n.up  <- n
                } else {
                    n.low <- n
                }
                # 'new' window width for the next iteration
                n <- ceiling( mean(c(n.low,n.up)) )    
                if( (n==n.up)|(n==n.low) ){
                    max.non.rejected <- max(which(H==0))
                    nt    <- all.widths.t[max.non.rejected]  # final window width
                    slope.t<- all.slopes.t[max.non.rejected]  # final slope estimate
                    level.t<- all.levels.t[max.non.rejected]  # final signal estimate
                    w <- c((length(y0)-nt+1):length(y0))
                    y <- y0[w]    # data vector for next step
                    M <- M0[w,w]  # matrix with slopes in next step
                    break
                }
            }
            w <- c((length(y0)-n+1):length(y0))
            y <- y0[w]    # data vector in next iteration
            M <- M0[w,w]  # matrix with slopes in next iteration
        } # end of width.search = "binary"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Set new window width to min.width
        if(width.search=="set.to.min"){
            if( n < 2*p.test || !reject.H0 ){
                nt   <- n
                break  # H0: "Fit is good" cannot be rejected or window width cannot be reduced
            } else {   # H0: "Fit is good" is rejected
                w <- (n-min.width+1):n
                n <- min.width     # -> set window width to minimum and re-estimate level and slope in smaller window
                y <- y[w]    # data vector in next iteration
                M <- M[w, w] # matrix with slopes in next iteration
            }
        } # end of linear one-step search



    } # end of repeat-loop for window width adaptation

return(list(y=y, level=level.t, slope=slope.t, width=nt, M=M,
            all.levels.t=all.levels.t, all.slopes.t=all.slopes.t, all.widths.t=all.widths.t, 
            I.t=I.t, nI=nI))

}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Internal Objects and Starting Values
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    n           <- start      # window width 
    I.win       <- 1:start    # indices in first window
    y           <- tseries[I.win]   # current data vector (observations in window)

    mu          <- double(N)  # all online RM levels iteratively estimated at each time t
    is.na(mu)   <- 1:N
    beta        <- double(N)  # all online RM slopes iteratively estimated at each time t
    is.na(beta) <- 1:N
    n.t         <- integer(N) # all window widths iteratively used at time each t
    is.na(n.t)  <- 1:start

    mu.list   <- vector("list", N) # list of levels
    beta.list <- vector("list", N) #         slopes
    n.t.list  <- vector("list", N) #         window widths at each time

    if(calc.qn) {
        qn          <- double(N) # vector of Qn scale estimates within each window
        is.na(qn)   <- 1:N
    }

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Matrix of RM slope within the first time window
    M            <- matrix(NA, ncol = n, nrow = n) # matrix with RM slopes
    for (k in 1:start){
        M[k, -k] <- (y[k] - y[-k]) / (I.win[k] - I.win[-k])
    }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# MAIN PROGRAMME LOOP OVER ALL TIME POINTS
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    for(t in start:N) {

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Evaluate adequate window width and Repeated Median estimates
        fit.n.t <- get.fit.and.n(y=y, M=M)

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Re-allocate internal objects
        y       <- fit.n.t$y
        n       <- fit.n.t$width
        nI      <- fit.n.t$nI
        M       <- fit.n.t$M
        I.t     <- fit.n.t$I.t
        
        # Save results
        n.t[t]  <- n
        beta[t] <- fit.n.t$slope

        # 'Restrict to range' rule for the estimated signal:
        if( !is.na(fit.n.t$level) & (rtr != 0) ) { 
            mu[t] <- restrict.mu(fit.n.t$level, y, I.t, rtr) 

        } else {
            mu[t]   <- fit.n.t$level
        }

        mu.list[[t]]   <- fit.n.t$all.levels.t
        beta.list[[t]] <- fit.n.t$all.slopes.t
        n.t.list[[t]]  <- fit.n.t$all.widths.t

        # Estimate scale (if wanted)
        if(calc.qn & (sum(!is.na(y[I.t])) >= ifelse( (nI < minNonNAs), nI, minNonNAs ))) {
            res <- get.res(y, fit.n.t$level, fit.n.t$slope, n)
            qn[t]   <- Qn(na.omit(res), finite.corr = TRUE)
        }

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Update step
        if(t != N) {
             if( ((n == max.width)&&is.numeric(max.width)) || sum(!is.na(y[I.t])) < ifelse( (nI < minNonNAs), nI, minNonNAs ) )  {
                # if window width maximal or not enough recent information
                # 'move' matrix of slope estimates and do not increase window width
                b.n     <- (tseries[t+1] - tseries[t+c((-n + 2):0)]) / (t+1 - (t+c((-n + 2):0)))
                M       <- cbind(rbind(M[-1, -1], b.n), c(b.n,NA))
                y       <- tseries[(t - n + 2):(t + 1)]
            } else {
                # enlarge matrix of slope estimates with new pairwise slopes and increase window width
                b.n     <- (tseries[t+1] - tseries[t+c((-n + 1):0)]) / (t+1 - (t+c((-n + 1):0)))
                M       <- cbind(rbind(M, b.n), c(b.n, NA))
                y       <- tseries[(t - n + 1):(t + 1)]
                n       <- n + 1
            }
        }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    } # end of MAIN-loop over all time points
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Extrapolation in the first time window 
    if(extrapolate){
    beta[1:(start-1)] <- beta[start]
    mu[1:(start-1)]   <- mu[start] - ((start-1):1)*beta[start]
    }

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Output
    result <- list(level=mu, slope=beta, width=n.t, 
                   level.list=mu.list, slope.list=beta.list, width.list=n.t.list, 
                   y=tseries, min.width=min.width, max.width=max.width, p.test=p.test,
                   minNonNAs=minNonNAs, calc.qn=calc.qn, rtr=rtr, extrapolate=extrapolate )

    if(calc.qn) {
        result$sigma      <- qn
    }
 
    return( structure( result , class="adore.filter") )
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default output
print.adore.filter <- function(x, ...) {
  N <- length( x$level )
  cat('$level \n')  
    if(N <= 100){
    print( x$level, ... )
  } else {
    print( x$level[1:(x$min.width + 4)] )
    cat('Only the first 5 online level estimations are printed.\n')    
  }
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default plot
plot.adore.filter <- function(x, ...) {
  if ( is.logical(x$max.width) && !x$max.width){
    Main <- paste("Online RM Filter with Adaptive Window Width\n(min width =", x$min.width, ", no max width)", sep="")
  } else {
    Main <- paste("Online RM Filter with Adaptive Window Width\n(min width=", x$min.width, ", max width=", x$max.width,")", sep="")
  }
  plot(x$y, type='l', main=Main, xlab='time', ylab='', ...)
  lines(x$level, col='red', lwd=2)
  legend(x="topleft", bty="n",legend=c("Time Series", "Filtered Signal"),lty=c(1, 1),col=c("black", "red"),lwd=c(1, 2))
}
