global.1m.ssc <- function(method, mean.diff, sd, cor, v = NULL, M = NULL, power = 0.8, alpha = 0.05) {

    if (missing(method)) stop("Missing 'method' argument.")
    if (missing(mean.diff)) stop("Missing 'mean.diff' argument.")
    if (missing(sd)) stop("Missing 'sd' argument.")
    if (missing(cor)) stop("Missing 'cor' argument.")	
    if (class(method) != "character") stop("The 'method' argument must be of type character.")
    if( (method != "Adj.Model") & (method != "Model") ) stop("The 'method' argument is misspecified.")	
    if ( (power < 0) | (power > 1)) stop("The 'power' argument should be between 0 and 1.")
    if ((alpha < 0) | (alpha > 1)) stop("The 'alpha' argument should be between 0 and 1.")
    if (!is.numeric(power)) stop("The 'power' argument should be numeric.")
    if (!is.numeric(alpha)) stop("The 'alpha' argument should be numeric.")      
    if (!is.vector(mean.diff)) stop("The 'mean.diff' argument should be a vector.")
    if (!is.vector(sd)) stop("The 'sd' argument should be a vector.")
    if (!is.matrix(cor)) stop("The 'cor' argument should be a matrix.")

        # Covariance matrix estimation
        Sigmahat <- (diag(sd) %*% cor %*% diag(sd))
    
    # Equation before (12), page 386 in our 2014 paper
    if (is.null(v) && is.null(M)) {
        coef <- 2
    } else {
        coef <- as.numeric(2 + (t(v) %*% solve(M) %*% v))
    }

    if (.f.n.global(10000, method = method, alpha = alpha, mean.diff = mean.diff, Sigmahat = Sigmahat, coef = coef, power = power) < 0) stop("The sample size is larger than 10,000.")

    # Research of the sample size
    n.value <- uniroot.integer(.f.n.global, interval = c(2, 10000), pos.side = TRUE, method = method, alpha = alpha, mean.diff = mean.diff, Sigmahat = Sigmahat, coef = coef, power = power)$root 
    
    return(n.value)
    
}

# Definition of power function
.powerfunc.global <- function(n, method, alpha, mean.diff, Sigmahat, coef) {
        
    m <- length(mean.diff)        

    if (method == "Adj.Model") {
        # Equation before (12), page 386 in our 2014 paper 
        W <- as.numeric((1 / n) * coef) * Sigmahat
    } else {
        # Line 2, page 396 in Appendix 2 in our 2014 paper
        W <- (2 / n) * Sigmahat
    }

    # Decentrality parameter. Equation (12) page 386 in our 2014 paper 
    ncp <- t(mean.diff) %*% solve(W) %*% mean.diff

    # Power 1 - \beta. Equation (14) page 387 in our 2014 paper:
    pow <- 1 - pchisq(qchisq(1 - alpha, df = m), df = m, ncp = ncp)

    return (pow)
}

.f.n.global <- function(n, method, alpha, mean.diff, Sigmahat, coef, power) {
    return(.powerfunc.global(n, method, alpha, mean.diff, Sigmahat, coef) - power)
}
