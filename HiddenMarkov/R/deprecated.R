"backward0.mmpp" <-
function(tau, Q, lambda){
    #   backward probabilities for MMPP
    #  tau contains the interevent times
    .Deprecated("forwardback.mmpp", package="HiddenMarkov",
          msg="'backward0.mmpp' is deprecated.
          Use 'forwardback.mmpp' instead, see help('forwardback.mmpp').")
    m <- nrow(Q)
    n <- length(tau)
    Lambda <- diag(lambda)
    logbeta <- matrix(rep(NA, m*(n+1)), nrow=(n+1))
    logbeta[(n+1),] <- 0
    phi <- matrix(rep(1/m, m), ncol=1)
    lscale <- log(m)
    decomp <- eigen(Q-Lambda, symmetric=FALSE)
    if (any(duplicated(decomp$values))) stop("repeated eigenvalues")
    S <- decomp$vectors
    Sinv <- solve(S)
    eigenval <- decomp$values
    for (i in seq(n, 1, -1)){
        phi <- S %*% diag(exp(eigenval*tau[i])) %*%
                 Sinv  %*% Lambda %*% phi
#       phi <- matrixexp((Q-Lambda)*tau[i]) %*% Lambda %*% phi
        logbeta[i,] <- log(phi) + lscale
        sumphi <- sum(phi)
        phi <- phi/sumphi
        lscale <- lscale + log(sumphi)
    }
    return(Re(logbeta))
}


"Baum.Welch" <-
function (x, Pi, delta, distn, pm, pn = NULL, nonstat = TRUE,
          maxiter = 500, tol = 1e-05, prt = TRUE,
          posdiff = (distn[1]!="glm")) 
{
    .Deprecated("BaumWelch", package="HiddenMarkov",
          msg="'Baum.Welch' is deprecated.
          Use 'BaumWelch' instead, see help('BaumWelch').")
    if (distn[1]!="glm"){
        Mstep <- parse(text=paste("Mstep.", distn,
                       "(x, cond, pm, pn)", sep=""))
    } else{
        Mstep <- parse(text=paste("Mstep.glm",
                  "(x, cond, pm, pn, distn[2], distn[3])", sep=""))
    }
    m <- nrow(Pi)
    n <- length(x)
    oldLL <- -Inf
    for (iter in 1:maxiter) {
        cond <- Estep(x, Pi, delta, distn, pm, pn)
        diff <- cond$LL - oldLL
        if (prt) {
            cat("iter =", iter, "\n")
            cat("LL =", formatC(cond$LL, digits=log10(1/tol)+2,
                                format="f"), "\n")
            cat("diff =", diff, "\n\n")
        }
        if (diff < 0 & posdiff) stop("Worse log-likelihood")
        if (abs(diff) < tol) break
        #----  Mstep  ----
        Pi <- diag(1/apply(cond$v, MARGIN = 2, FUN = sum)) %*% 
            apply(cond$v, MARGIN = c(2, 3), FUN = sum)
        if (nonstat) delta <- cond$u[1, ]
        else delta <- compdelta(Pi)
        pm <- eval(Mstep)
        oldLL <- cond$LL
    }
    return(list(delta = delta, Pi = Pi, u = cond$u, 
           v = cond$v, pm=pm, LL = cond$LL,
           iter = iter, diff = diff))
}


"Baum.Welch0.mmpp" <-
function (tau, Q, delta, lambda, nonstat = TRUE,
          maxiter = 500, tol = 1e-05, prt = TRUE,
          converge=expression(diff < tol)) 
{
    #  Using the method in Ryden (1996) - unscaled
    .Deprecated("BaumWelch", package="HiddenMarkov",
          msg="'Baum.Welch0.mmpp' is deprecated.
          Use 'BaumWelch' instead, see help('BaumWelch').")
    m <- nrow(Q)
    n <- length(tau)
    oldLL <- -Inf
    for (iter in 1:maxiter) {
        cond <- Estep0.mmpp(tau, Q, delta, lambda)
#       overflow problems
#       LL <- log(sum(exp(cond$logalpha[(n+1),])))
        LL <- logLikmmpp(tau, Q, delta, lambda)
        diff <- LL - oldLL
        if (prt) {
            cat("iter =", iter, "\n")
            cat("LL =", formatC(LL, digits=log10(1/tol)+2,
                                format="f"), "\n")
            cat("diff =", diff, "\n\n")
        }
        if (diff < 0){
            stop("Worse log-likelihood on last iteration")
        }
        if (eval(converge)) break
        #----  Mstep  ----
        Q <- Q * (diag(1/diag(cond$A)) %*% cond$A)
        diag(Q) <- 0
        diag(Q) <- -apply(Q, MARGIN=1, FUN=sum)
        lambda <- cond$B/diag(cond$A)
        if (nonstat) delta <- exp(cond$logalpha[1, ] +
                                  cond$logbeta[1, ] - LL)
        else delta <- compdelta(solve(diag(lambda) -
                                Q) %*% diag(lambda))
        oldLL <- LL
    }
    return(list(delta = delta, Q = Q, lambda = lambda, 
           LL = LL, iter = iter, diff = diff))
}


"Baum.Welch.mmpp" <-
function (tau, Q, delta, lambda, nonstat = TRUE,
          maxiter = 500, tol = 1e-05, prt = TRUE,
          converge=expression(diff < tol)) 
{
    #   scaled version of Ryden (1996)
    .Deprecated("BaumWelch", package="HiddenMarkov",
          msg="'Baum.Welch.mmpp' is deprecated.
          Use 'BaumWelch' instead, see help('BaumWelch').")
    m <- nrow(Q)
    n <- length(tau)
    oldLL <- -Inf
    for (iter in 1:maxiter) {
        cond <- Estep.mmpp(tau, Q, delta, lambda)
        LL <- cond$LL
        diff <- LL - oldLL
        if (prt) {
            cat("iter =", iter, "\n")
            cat("LL =", formatC(LL, digits=log10(1/tol)+2,
                                format="f"), "\n")
            cat("diff =", diff, "\n\n")
        }
        if (diff < 0){
            stop("Worse log-likelihood on last iteration")
        }
        if (eval(converge)) break
        #----  Mstep  ----
        Q <- Q * (diag(1/diag(cond$A)) %*% cond$A)
        diag(Q) <- 0
        diag(Q) <- -apply(Q, MARGIN=1, FUN=sum)
        lambda <- cond$B/diag(cond$A)
        if (nonstat) delta <- exp(cond$logalpha[1, ] +
                                  cond$logbeta[1, ])
        else delta <- compdelta(solve(diag(lambda) -
                                Q) %*% diag(lambda))
        oldLL <- LL
    }
    return(list(delta = delta, Q = Q, lambda = lambda,
           LL = LL, iter = iter, diff = diff))
}


"forward0.mmpp" <-
function(tau, Q, delta, lambda){
    #  forward probs for MMPP
    #  tau contains the interevent times
    .Deprecated("forwardback.mmpp", package="HiddenMarkov",
          msg="'forward0.mmpp' is deprecated.
          Use 'forwardback.mmpp' instead, see help('forwardback.mmpp').")
    m <- nrow(Q)
    n <- length(tau)
    Lambda <- diag(lambda)
    phi <- matrix(delta, nrow=1)
    logalpha <- matrix(rep(NA, m*(n+1)), nrow=(n+1))
    logalpha[1,] <- log(phi)
    lscale <- 0
    decomp <- eigen(Q-Lambda, symmetric=FALSE)
    if (any(duplicated(decomp$values))) stop("repeated eigenvalues")
    S <- decomp$vectors
    Sinv <- solve(S)
    eigenval <- decomp$values
    for (i in 2:(n+1)){
        phi <- phi %*% S %*% diag(exp(eigenval*tau[i-1])) %*%
                 Sinv %*% Lambda
#       phi <- phi %*% matrixexp((Q-Lambda)*tau[i-1]) %*% Lambda
        sumphi <- sum(phi)
        phi <- phi/sumphi
        lscale <- lscale + log(sumphi)
        logalpha[i,] <- log(phi) + lscale
    }
    return(Re(logalpha))
}


"sim.hmm" <-
function (n, initial, Pi, distn, pm, pn = NULL) 
{
    .Deprecated("simulate", package="HiddenMarkov",
          msg="'sim.hmm' is deprecated.
          Use 'simulate' instead, see help('simulate').")
    rname <- paste("r", distn, sep="")
    y <- sim.markov(n, initial, Pi)
    x <- rep(NA, n)
    for (i in 1:n) {
        x[i] <- do.call(rname, c(list(n=1), getj(pm, y[i]),
                        getj(pn, i)))
    }
    return(list(x = x, y = y))
}


"sim.hmm1" <-
function (n, initial, Pi, distn, pm) 
{
    .Deprecated("simulate", package="HiddenMarkov",
          msg="'sim.hmm1' is deprecated.
          Use 'simulate' instead, see help('simulate').")
    rfunc <- eval(parse(text=paste("r", distn, sep="")))
    nms <- names(pm)
    newargs <- paste(nms, "=", "pm[[", seq(1,length(nms)),
                     "]][k]", sep="", collapse=",")
    newargs <- paste("alist(n=, ", newargs, ")", sep="")
    formals(rfunc) <- eval(parse(text=newargs))
    y <- sim.markov(n, initial, Pi)
    x <- rep(NA, n)
    m <- nrow(Pi)
    for (k in 1:m) {
        a <- (y == k)
        x[a] <- rfunc(n=sum(a))
    }
    return(list(x = x, y = y))
}


"sim.markov" <-
function(n, initial, Pi){
    .Deprecated("simulate", package="HiddenMarkov",
          msg="'sim.markov' is deprecated.
          Use 'simulate' instead, see help('simulate').")
    #    simulate a Markov Chain
    x <- rep(NA, n)
    x[1] <- initial
    m <- nrow(Pi)
    for (i in 2:n)
        x[i] <- sample(x=1:m, size=1, prob=Pi[(x[i-1]),])
    return(x)
}


"sim.mmpp" <-
function (n, initial, Q, lambda) 
{
    #   ********   VALID FOR >=2 STATES   ********
    .Deprecated("BaumWelch", package="HiddenMarkov",
          msg="'sim.mmpp' is deprecated.
          Use 'simulate' instead, see help('simulate').")
    #   y    contains sequence of Markov states
    #   x    transition time to next state
    #   tau  times of Poisson events
    m <- ncol(Q)
    Pi <- diag(m) - diag(1/diag(Q)) %*% Q
    ys <- rep(NA, n+1)
    tau <- rep(NA, n+1)
    #    the length of x and y may be too short
    #    gets extended later if required
    x <- rep(NA, n*10)
    y <- rep(NA, n*10)
    y[1] <- ys[1] <- initial
    x[1] <- tau[1] <- 0
    i <- j <- 2
    while (TRUE){
        #   sim time spent in Markov state y[i-1]
        y[i] <- sample(x=1:m, size=1, prob=Pi[(y[i-1]),])
        x[i] <- x[i-1] + rexp(1, rate=-Q[y[i-1], y[i-1]])
        t0 <- x[i-1]
        while(TRUE){
            #   sim times of Poisson events
            ti <- t0 + rexp(1, rate=lambda[y[i-1]])
            if (ti < x[i]){
                tau[j] <- t0 <- ti
                ys[j] <- y[i-1]
                j <- j + 1
                if (j==n+2)
                    return(list(x=x[1:i], y=y[1:i], tau=tau, ys=ys))
            }
            else break
        }
        i <- i+1
        #    extend x and y if too short
        if (i > length(x)){
            x <- c(x, rep(NA, n*10))
            y <- c(y, rep(NA, n*10))
        }
    }
}


"Viterbihmm" <-
function (x, Pi, delta, distn, pm, pn = NULL) 
{
    .Deprecated("Viterbi", package="HiddenMarkov",
          msg="'Viterbihmm' is deprecated.
          Use the dthmm model object with 'Viterbi', see help('Viterbi').")
    dfunc <- makedensity(distn)
    n <- length(x)
    m <- nrow(Pi)
    nu <- matrix(NA, nrow = n, ncol = m)
    y <- rep(NA, n)
    nu[1, ] <- log(delta) + dfunc(x=x[1], pm, getj(pn, 1),
                                  log=TRUE)
    logPi <- log(Pi)
    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        nu[i, ] <- apply(matrixnu + logPi, 2, max) +
                      dfunc(x=x[i], pm, getj(pn, i),
                                  log=TRUE)
    }
    if (any(nu[n, ] == -Inf)) 
        stop("Problems With Underflow")
    y[n] <- which.max(nu[n, ])
    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[, y[i + 
        1]] + nu[i, ])
    return(y)
}

"logLikmmpp" <-
function(tau, Q, delta, lambda){
    #   Log-likelihood of MMPP
    .Deprecated("logLik", package="HiddenMarkov",
          msg="'logLikmmpp' is deprecated.
          Use 'logLik' instead, see help('logLik').")
    n <- length(tau)
    Lambda <- diag(lambda)
    phi <- delta
    LL <- 0
    decomp <- eigen(Q-Lambda, symmetric=FALSE)
    if (any(duplicated(decomp$values))) stop("repeated eigenvalues")
    S <- decomp$vectors
    post <- solve(S) %*% Lambda
    eigenval <- decomp$values
    for (i in 1:n){
        phi <- phi %*% S %*% diag(exp(eigenval*tau[i])) %*%
                 post
        sumphi <- sum(phi)
        LL <- LL + log(sumphi)
        phi <- phi/sumphi
    }
    return(LL)
}


