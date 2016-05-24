"dat2bernqua" <-
function(f, x, bern.control=NULL,
               poly.type=c("Bernstein", "Kantorovich", "Cheng", "Parzen"),
               bound.type=c("none", "sd", "Carv", "either"),
               fix.lower=NULL, fix.upper=NULL, p=0.05, listem=FALSE) {

    poly.type <- match.arg(poly.type)
    if(poly.type == "Parzen") {
       f <- 0.5 # just an arbitrary value to bypass the check.fs
       # The Parzen method uses its own specific definition of the f values
    }


    if(! is.null(bern.control)) {
         poly.type  <- bern.control$poly.type
         bound.type <- bern.control$bound.type
         fix.lower  <- bern.control$fix.lower
         fix.upper  <- bern.control$fix.upper
         p          <- bern.control$p
    }

    if(p < 1E-6 || p >= (1 - 1E-6)) {
        warning("p is too small or too large (ad hoc decision, see source code), returning NA")
        return(NA)
    }
    poly.type  <- match.arg(poly.type)
    bound.type <- match.arg(bound.type)
    if(! check.fs(f)) return()
    x <- sort(x); n <- length(x)
    if(! is.null(fix.lower) && x[1] < fix.lower) {
       warning("The observed minimum is less than the declared lower bounds, resetting to observed minimum")
       fix.lower <- x[1]
    }
    if(! is.null(fix.upper) && x[n] > fix.upper) {
       warning("The observed maximum is greater than the declared upper bounds, resetting to observed maximum")
       fix.upper <- x[n]
    }
    lam2 <- lmoms(x, nmom=2)$lambdas[2]
    # Compute sd-based bounds
    sd.lower <- x[1] - lam2*sqrt(pi/n)
    sd.upper <- x[n] + lam2*sqrt(pi/n)
    # Compute de Carvalho bounds
    a <- (1-p)^(-2) - 1
    Carv.lower <- x[1] - (x[2] - x[1])/a
    Carv.upper <- x[n] + (x[n] - x[(n-1)])/a

    if(poly.type == "Parzen") {
       fix.lower <- ifelse(is.null(fix.lower), x[1], fix.lower)
       fix.lower <- min(fix.lower, x[1])
       x[n+1] <- fix.lower; x <- sort(x)
       u <- sapply(1:n, function(r) { return((r-1)/n) })
       qua <- sapply(1:n, function(r) { A <- (n*((r/n) - u[r]))*x[(r-1)+1]; B <- n*(u[r] - (r-1)/n)*x[(r+1)]; return(A+B) }) 
       return(list(f=u, x=qua))
    }

    if(bound.type == "sd") {
        #message("Using the standard deviation support if either end is larger (smaller) than the respective data minimum (maximum)")
        fix.lower <- max(c(sd.lower, fix.lower))
        fix.upper <- min(c(sd.upper, fix.upper))
    } else if(bound.type == "Carv") {
        #message("Using the de Carvalho support if either end is larger (smaller) than the respective data minimum (maximum)")
        fix.lower <- max(c(Carv.lower, fix.lower))
        fix.upper <- min(c(Carv.upper, fix.upper))
    } else if(bound.type == "either") {
        #message("Using either the standard deviation or de Carvalho support if either end is larger (smaller) than the respective data minimum (maximum)")
        fix.lower <- max(c(sd.lower, Carv.lower, fix.lower))
        fix.upper <- min(c(sd.upper, Carv.upper, fix.upper))
    } else if(bound.type == "none") {
        if(is.null(fix.lower)) fix.lower <- x[1]
        if(is.null(fix.upper)) fix.upper <- x[n]
    } else {
       # Do nothing
    }

    qua <- vector(mode="numeric", length=length(f))
    for(i in 1:length(f)) {
       myf  <- log(f[i]);     if(! is.finite(myf))  myf  <- 0
       myfc <- log(1 - f[i]); if(! is.finite(myfc)) myfc <- 0
       if(poly.type == "Cheng") {
             tmp <- sapply(1:n, function(k) {
                        xk <- x[k]
                        return(xk * exp(lchoose(n-1,k-1) + (k-1)*myf + (n-k)*myfc)) })
             qua[i] <- sum(tmp)
       } else {
          if(poly.type == "Bernstein") {
             tmp <- sapply(0:(n+1), function(k) {
                        xk <- x[k]
                        if(k ==     0) xk <- fix.lower
                        if(k == (n+1)) xk <- fix.upper
                        return(xk * exp(lchoose(n+1,k) + k*myf + (n+1-k)*myfc)) })
             qua[i] <- sum(tmp)
          } else if(poly.type == "Kantorovich") {
             tmp <- sapply(0:n, function(k) {
                        xk   <- ifelse(k == 0, fix.lower, x[k]  )
                        xkp1 <- ifelse(k == n, fix.upper, x[k+1])
                        return((xk+xkp1) * exp(lchoose(n,k) + k*myf + (n-k)*myfc)) })
             qua[i] <- sum(tmp)/2
          } else {
             stop("Should not be here in the logic")
          }
       }
    }

    if(listem) {
       z <- list(f=f, x=qua, n=n, p=p, fix.lower=fix.lower, fix.upper=fix.upper)
       return(z)
    } else {
       return(qua)
    }
}


