# Nelder-Mead with Box constraints
nmkb2 <- function (par, fn, lower = -Inf, upper = Inf, control = list(), ...) 
{
  #DEBUG=FALSE
    ctrl <- list(tol = 1e-06, maxfeval = min(5000, max(1500, 
        20 * length(par)^2)), regsimp = TRUE, maximize = FALSE, 
        restarts.max = 3, trace = FALSE)
    namc <- match.arg(names(control), choices = names(ctrl), 
        several.ok = TRUE)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    if (!is.null(names(control))) 
        ctrl[namc] <- control
    ftol <- ctrl$tol
    maxfeval <- ctrl$maxfeval
    regsimp <- ctrl$regsimp
    restarts.max <- ctrl$restarts.max
    maximize <- ctrl$maximize
    trace <- ctrl$trace
    n <- length(par)
  
  
  ##
  ## IMPORTANT NOTE: Nelder-Mead crashes, if the starting value is right at the boundary
  ## in any dimension (i.e. the distance to the boundary is < 1e-14). 
  ## Since the optimum for some problems is precisely at the boundary, 
  ## this is a problem. We cure this symptomatically by shifting the starting point 
  ## a little bit (1e-13) 'inside' in any dimension.
  ##
  clipToRegion <- function(xc,lowerP,upperP) {
    xc = pmax(xc,lowerP+1e-13)
    xc = pmin(xc,upperP-1e-13)
    return(xc)
  }
  par = clipToRegion(par,lower,upper)  

g <- function(x) {
gx <- x
gx[c1] <- atanh(2 * (x[c1] - lower[c1]) / (upper[c1] - lower[c1]) - 1)   
gx[c3] <- log(x[c3] - lower[c3])
gx[c4] <- log(upper[c4] - x[c4])
gx
}

ginv <- function(x) {
gix <- x
gix[c1] <- lower[c1] + (upper[c1] - lower[c1])/2 * (1 + tanh(x[c1]))
gix[c3] <- lower[c3] + exp(x[c3])
gix[c4] <- upper[c4] - exp(x[c4])
gix
}

if (length(lower) == 1) lower <- rep(lower, n)
if (length(upper) == 1) upper <- rep(upper, n)

if (any(c(par < lower, upper < par))) stop("Infeasible starting values!", call.=FALSE)

low.finite <- is.finite(lower)
upp.finite <- is.finite(upper)
c1 <- low.finite & upp.finite  # both lower and upper bounds are finite 
c2 <- !(low.finite | upp.finite) # both lower and upper bounds are infinite
c3 <- !(c1 | c2) & low.finite # finite lower bound, but infinite upper bound
c4 <- !(c1 | c2) & upp.finite  # finite upper bound, but infinite lower bound
    
if (all(c2)) stop("Use `nmk()' for unconstrained optimization!", call.=FALSE)

    if (maximize) 
        fnmb <- function(par) -fn(ginv(par), ...)
    else fnmb <- function(par) fn(ginv(par), ...)
    
    x0 <- g(par)
    if (n == 1) 
        stop(call. = FALSE, "Use `optimize' for univariate optimization")
    if (n > 30) 
        warning("Nelder-Mead should not be used for high-dimensional optimization")
    V <- cbind(rep(0, n), diag(n))
    f <- rep(0, n + 1)
    f[1] <- fnmb(x0)
    V[, 1] <- x0
    scale <- max(1, sqrt(sum(x0^2)))
    if (regsimp) {
      #browser()
      alpha <- scale/(n * sqrt(2)) * c(sqrt(n + 1) + n - 1, 
            sqrt(n + 1) - 1)
        V[, -1] <- (x0 + alpha[2])
        diag(V[, -1]) <- x0[1:n] + alpha[1]
        for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
    }
    else {
      #browser()
      V[, -1] <- x0 + scale * V[, -1]
        for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
    }
    #if (DEBUG && any(is.na(f))) browser()
    f[is.nan(f)] <- Inf
    f[is.na(f)] <- Inf    # /WK/
    nf <- n + 1
    ord <- order(f)
    f <- f[ord]
    V <- V[, ord]
    rho <- 1
    gamma <- 0.5
    chi <- 2
    sigma <- 0.5
    conv <- 1
    oshrink <- 0
    restarts <- 0
    orth <- 0
    dist <- f[n + 1] - f[1]
    v <- V[, -1] - V[, 1]
    delf <- f[-1] - f[1]
    diam <- sqrt(colSums(v^2))
    sgrad <- c(crossprod(t(v), delf))
    alpha <- 1e-04 * max(diam)/sqrt(sum(sgrad^2))
    simplex.size <- sum(abs(V[, -1] - V[, 1]))/max(1, sum(abs(V[, 
        1])))
    #if (DEBUG && any(is.na(simplex.size))) browser()
    if (is.nan(simplex.size)) simplex.size=Inf    # /WK/ It can happen that the initial V
    if (is.na(simplex.size)) simplex.size=Inf     # produces an (Inf/Inf = NaN)-situation. 
                                                  # Then the while-condition below will crash.
                                                  # We cure this symptomatically with 
                                                  #     simplex.size=Inf
    #if (DEBUG && any(is.na(dist))) browser()
    if (is.nan(dist)) dist <- Inf     # /WK/
    if (is.na(dist)) dist <- Inf      # /WK/
    itc <- 0
    conv <- 0
    message <- "Succesful convergence"
    #browser()
    while (nf < maxfeval & restarts < restarts.max & dist > ftol & 
        simplex.size > 1e-06) {
        fbc <- mean(f)
        happy <- 0
        itc <- itc + 1
        xbar <- rowMeans(V[, 1:n])
        xr <- (1 + rho) * xbar - rho * V[, n + 1]
        fr <- fnmb(xr)
        nf <- nf + 1
        #if (DEBUG && any(is.na(fr))) browser()
        if (is.nan(fr))  fr <- Inf
        if (is.na(fr))  fr <- Inf   # /WK/
        if (fr >= f[1] & fr < f[n]) {
            happy <- 1
            xnew <- xr
            fnew <- fr
        }
        else if (fr < f[1]) {
            xe <- (1 + rho * chi) * xbar - rho * chi * V[, n + 
                1]
            fe <- fnmb(xe)
            #if (DEBUG && any(is.na(fe))) browser()
            if (is.nan(fe))   fe <- Inf
            if (is.na(fe))   fe <- Inf    # /WK/
            nf <- nf + 1
            if (fe < fr) {
                xnew <- xe
                fnew <- fe
                happy <- 1
            }
            else {
                xnew <- xr
                fnew <- fr
                happy <- 1
            }
        }
        else if (fr >= f[n] & fr < f[n + 1]) {
            xc <- (1 + rho * gamma) * xbar - rho * gamma * V[, 
                n + 1]
            fc <- fnmb(xc)
            #if (DEBUG && any(is.na(fc))) browser()    
            if (is.nan(fc))  fc <- Inf
            if (is.na(fc))  fc <- Inf     # /WK/
            nf <- nf + 1
            if (fc <= fr) {
                xnew <- xc
                fnew <- fc
                happy <- 1
            }
        }
        else if (fr >= f[n + 1]) {
            xc <- (1 - gamma) * xbar + gamma * V[, n + 1]
            fc <- fnmb(xc)
            #if (DEBUG && any(is.na(fc))) browser()
            
            if (is.nan(fc))  fc <- Inf
            if (is.na(fc))  fc <- Inf     # /WK/
            nf <- nf + 1
            if (fc < f[n + 1]) {
                xnew <- xc
                fnew <- fc
                happy <- 1
            }
        }
        if (happy == 1 & oshrink == 1) {
            fbt <- mean(c(f[1:n], fnew))
            delfb <- fbt - fbc
            armtst <- alpha * sum(sgrad^2)
            if (delfb > -armtst/n) {
                if (trace) 
                  cat("Trouble - restarting: \n")
                restarts <- restarts + 1
                orth <- 1
                diams <- min(diam)
                sx <- sign(0.5 * sign(sgrad))
                happy <- 0
                V[, -1] <- V[, 1]
                diag(V[, -1]) <- diag(V[, -1]) - diams * sx[1:n]
            }
        }
        if (happy == 1) {
            V[, n + 1] <- xnew
            f[n + 1] <- fnew
            ord <- order(f)
            V <- V[, ord]
            f <- f[ord]
        }
        else if (happy == 0 & restarts < restarts.max) {
            if (orth == 0) 
                orth <- 1
            V[, -1] <- V[, 1] - sigma * (V[, -1] - V[, 1])
            for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
            nf <- nf + n
            ord <- order(f)
            V <- V[, ord]
            f <- f[ord]
        }
        v <- V[, -1] - V[, 1]
        delf <- f[-1] - f[1]
        diam <- sqrt(colSums(v^2))
        simplex.size <- sum(abs(v))/max(1, sum(abs(V[, 1])))
        if (is.nan(simplex.size)) simplex.size=Inf    # /WK/ It can happen that V
        if (is.na(simplex.size)) simplex.size=Inf     # produces an (Inf/Inf = NaN)-situation. 
                                                      # Then the while-condition above will crash.
                                                      # We cure this symptomatically with 
                                                      #     simplex.size=Inf
        #if (DEBUG && any(is.na(f))) browser()
        
        f[is.nan(f)] <- Inf
        f[is.na(f)] <- Inf    # /WK/
        dist <- f[n + 1] - f[1]
        sgrad <- c(crossprod(t(v), delf))
        if (trace & !(itc%%2)) 
            cat("iter: ", itc, "\n", "value: ", f[1], "\n")
        
        #if (DEBUG && any(is.na(dist))) browser()        
        if (is.nan(dist)) dist <- Inf     # /WK/
        if (is.na(dist)) dist <- Inf      # /WK/
        #browser()
    }
    if (dist <= ftol | simplex.size <= 1e-06) {
        conv <- 0
        message <- "Successful convergence"
    }
    else if (nf >= maxfeval) {
        conv <- 1
        message <- "Maximum number of fevals exceeded"
    }
    else if (restarts >= restarts.max) {
        conv <- 2
        message <- "Stagnation in Nelder-Mead"
    }
    return(list(par = ginv(V[, 1]), value = f[1] * (-1)^maximize, feval = nf, 
        restarts = restarts, convergence = conv, message = message))
}
