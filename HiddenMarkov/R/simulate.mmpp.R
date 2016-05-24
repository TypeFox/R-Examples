simulate.mmpp <- function(object, nsim=1, seed=NULL, ...){
    #   VALID FOR >=2 STATES
    #   y    contains sequence of Markov states
    #   x    transition time to next state
    #   tau  times of Poisson events
    if (!is.null(seed)) set.seed(seed)
    m <- ncol(object$Q)
    #------------------------
    if (sum(object$delta)!=1) stop("Invalid delta")
    if (any(object$delta==1))
        initial <- (1:m)[as.logical(object$delta)]
    else
        initial <- sample(m, 1, prob=object$delta)
    #------------------------
    Q <- object$Q
    lambda <- object$lambda
    m <- ncol(Q)
    Pi <- diag(m) - diag(1/diag(Q)) %*% Q
    ys <- rep(NA, nsim+1)
    tau <- rep(NA, nsim+1)
    #    the length of x and y may be too short
    #    gets extended later if required
    x <- rep(NA, nsim*10)
    y <- rep(NA, nsim*10)
    y[1] <- ys[1] <- initial
    x[1] <- tau[1] <- 0
    i <- 1
    j <- 2
    while (j < nsim+2){
        i <- i+1
        #    extend x and y if too short
        if (i > length(x)){
            x <- c(x, rep(NA, nsim*10))
            y <- c(y, rep(NA, nsim*10))
        }
        #   sim time spent in Markov state y[i-1]
        y[i] <- sample(x=1:m, size=1, prob=Pi[(y[i-1]),])
        x[i] <- x[i-1] + rexp(1, rate=-Q[y[i-1], y[i-1]])
        t0 <- x[i-1]
        while(j < nsim+2){
            #   sim times of Poisson events
            ti <- t0 + rexp(1, rate=lambda[y[i-1]])
            if (ti < x[i]){
                tau[j] <- t0 <- ti
                ys[j] <- y[i-1]
                j <- j + 1
            }
            else break
        }
    }
    object$x <- x[1:i]
    object$y <- y[1:i]
    object$tau <- tau
    object$ys <- ys
    return(object)
}

