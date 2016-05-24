# This file includes following functions
# cm.matrix, cm.quantile, cm.cs, cm.ref, cm.rnorm, cm.rnorm.cor, cm.state,
# cm.val, cm.portfolio, cm.gain, cm.CVaR, cm.hist


# Function for testing if M is a migration matrix

"cm.matrix" <- function(M)
{
    if (!is.matrix(M)) 
        stop("M is not a matrix")
    if (dim(M)[1] < 2 || dim(M)[2] < 2 || dim(M)[1] != dim(M)[2]) 
        stop("incorrect dimension of migration matrix")
    if (length(M[M < 0]) != 0) 
        stop("negative value in input data") 
    if (length(M[M > 1]) != 0) 
        stop("probabilities larger than 1 in input data")
    
    # Test if the sum of each row is approximatively 1
    y <- numeric(dim(M)[1])
    for (i in 1:dim(M)[1]) {
        y[i] <- sum(M[i,])
    }

    y <- round(y, digits = 12)
    
    if (length(y[y != 1]) != 0)
        stop("the sum of each row is not equal 1")
}


# Computation of thresholds

"cm.quantile" <- function(M)
{
    cm.matrix(M)
    mpRev <- M[1:(dim(M)[1] - 1), dim(M)[2]:1]
    cumMPRev <- t(apply(mpRev, 1, cumsum))
    q <- qnorm(cumMPRev)
    
    # Adding of negative threshold
    q <- cbind(-Inf, q)

    return(q)
}


# Computation of credit spreads

"cm.cs" <- function(M, lgd)
{
    cm.matrix(M)
    
    if (lgd < 0 || lgd > 1)
        stop("loss given default not between 0 and 1")
    pd <- M[1:(dim(M)[1]-1), dim(M)[2]]
    cs <- -(log(1 - lgd * pd))
    
    return(cs)
}


# Computation of reference value

"cm.ref" <- function(M, lgd, ead, r, rating)
{
    cm.matrix(M)
    
    if (r < 0 || r > 1)
        stop("interest rate not between 0 and 1") 
    if (length(ead) != length(rating))
        stop("vectors ead and rating should have the same length")
    cs <- cm.cs(M, lgd)
    
    # Value of each position without rating change
    y.constVal <- ead * exp(-(r + cs[rating]))
    
    # Portfolio value without rating change at time t = 1
    y.constPV <- sum(y.constVal)

    list(constVal = y.constVal, constPV = y.constPV) 
}


# Simulation of random numbers with antithetic sampling

"cm.rnorm" <- function(N, n)
{
    if (N <= 0 || length(N) != 1)
        stop("N should be greater 0 and its length should be 1")
    if ((N <=1 && n <=1) || N > n)
        stop("N and n have no valid input data")
    X <- matrix(rnorm(N * n / 2), N, n / 2)
    Y <- cbind(X, -X)

    return(Y)   
}


# Simulation of correlated random numbers

"cm.rnorm.cor" <- function(N, n, rho)
{
    if (!is.matrix(rho)) 
        stop("rho is not a matrix")
        
    y <- eigen(rho)$values
    if (length(y[y < 0]) != 0)
        stop("rho should be positiv definit or positiv semidefinit")
        
    # Cholesky matrix
    A <- t(chol(rho))
    Y <- cm.rnorm(N, n)
    Y <- A %*% Y

    return(Y)
}


# Computation of state space: the state space is computed at time t = 1, 
# that is credit positions of each company and for each migration

"cm.state" <- function(M, lgd, ead, N, r)
{
    cm.matrix(M)
    
    if (r < 0 || r > 1)
        stop("interest rate not between 0 and 1") 
    if (length(ead) != N)
        stop("the length of vector ead and N should be equal")
        
    cs <- cm.cs(M, lgd)
    cs2 <- matrix(rep(cs, N), N, dim(M)[1] - 1, byrow=T)
    ead2 <- matrix(rep(ead, dim(M)[1] - 1), N, dim(M)[1] - 1, byrow=F)
    V <- ead2 * exp(-(r + cs2))
    
    # Special case: event of default
    V <- cbind(V, ead * (1 - lgd))
    
    return(V)
}


# Valuation of credit positions of each scenario

"cm.val" <- function(M, lgd, ead, N, n, r, rho, rating)
{
    cm.matrix(M)
    
    if (length(ead) != length(rating))
        stop("vectors ead and rating should have the same length")
    if (!is.matrix(rho)) 
        stop("rho is not a matrix")
        
    y <- eigen(rho)$values
    if (length(y[y < 0]) != 0)
        stop("rho should be positiv definit or positiv semidefinit")
    
    # Allocation in rating classes and identification of each credit position
    V <- cm.state(M, lgd, ead, N, r)
    simV <- matrix(0,N,n)
    q <- cm.quantile(M)
    Y <- cm.rnorm.cor(N, n, rho)
    for (i in 1:N) {
        l <- q[rating[i], ]
        simClasses <- findInterval(Y[i, ], l)
        simClasses <- (simClasses - (dim(M)[1] + 1)) * (-1)
        simV[i,] <- V[i,simClasses]
    }

    return(simV) 
}


# Computation of simulated portfolio values

"cm.portfolio" <- function(M, lgd, ead, N, n, r, rho, rating)
{
    cm.matrix(M)
    
    simV <- cm.val(M, lgd, ead, N, n, r, rho, rating)
    simPV <- colSums(simV)

    return(simPV)
}


# Computation of simulated profits and losses

"cm.gain" <- function(M, lgd, ead, N, n, r, rho, rating)
{
    cm.matrix(M)
    
    constPV <- cm.ref(M, lgd, ead, r, rating)$constPV
    simPV <- cm.portfolio(M, lgd, ead, N, n, r, rho, rating)

    simGV <- simPV - constPV

    return(simGV)
}


# Computation of credit value at risk

"cm.CVaR" <- function(M, lgd, ead, N, n, r, rho, alpha, rating)
{
    cm.matrix(M)
    
    if (alpha < 0 || alpha > 1)
        stop("confidence level alpha not between 0 and 1")
    
    simGV <- cm.gain(M, lgd, ead, N, n, r, rho, rating)
    CVaR <- -quantile(simGV, 1 - alpha)

    return(CVaR)
}


# Profit / Loss Distribution as histogram

"cm.hist" <- function(M, lgd, ead, N, n, r, rho, rating, 
                      col = "steelblue4", main = "Profit / Loss Distribution", 
                      xlab = "profit / loss", ylab = "frequency")
{
    cm.matrix(M)
    
    simGV <- cm.gain(M, lgd, ead, N, n, r, rho, rating)
    b <- seq(min(simGV), max(simGV), length=(max(simGV)-min(simGV))/(2*n))
    
    hist(simGV, main = main, xlab = xlab, ylab = ylab, breaks = b, col = col)
}