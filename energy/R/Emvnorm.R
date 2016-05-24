mvnorm.etest <- 
function(x, R = 999) 
{
    # parametric bootstrap E-test for multivariate normality
    if (is.vector(x)) {
        n <- length(x)
        d <- 1
        bootobj <- boot(x, statistic = normal.e, R = R, sim = "parametric", 
            ran.gen = function(x, y) {return(rnorm(n)) })
        }
        else {
        n <- nrow(x)
        d <- ncol(x)
        bootobj <- boot(x, statistic = mvnorm.e, R = R, sim = "parametric", 
            ran.gen = function(x, y) {
                return(matrix(rnorm(n * d), nrow = n, ncol = d)) })
        }
    p <- 1 - mean(bootobj$t < bootobj$t0)
    names(bootobj$t0) <- "E-statistic"
    e <- list(statistic = bootobj$t0,
              p.value = p,
              method = "Energy test of multivariate normality: estimated parameters",
              data.name = paste("x, sample size ", n, ", dimension ", d, ", replicates ", R, sep = ""))
    class(e) <- "htest"        
    e                 
}

mvnorm.e <- 
function(x) 
{
    # E-statistic for multivariate normality
    if (is.vector(x)) return(normal.e(x))
    n <- nrow(x)
    d <- ncol(x)
    if (n < 2) return(normal.e(x))
    z <- scale(x, scale = FALSE)    #subtract column means and 
    ev <- eigen(var(x), symmetric = TRUE)    #compute S^(-1/2)
    P <- ev$vectors
    lambda <- ev$values    
    y <- z %*% (P %*% diag(1 / sqrt(lambda)) %*% t(P))
    if (any(!is.finite(y))) return (NA)
    stat <- 0
    e <- .C("mvnEstat", y = as.double(t(y)), byrow = as.integer(TRUE),
            nobs = as.integer(n), dim = as.integer(d), 
            stat = as.double(stat), PACKAGE = "energy")$stat
    e
}

normal.e <- 
function(x) 
{
   x <- as.vector(x)
   y <- sort(x)
   n <- length(y)
   if (y[1] == y[n]) return (NA)
   y <- scale(y) 
   K <- seq(1 - n, n - 1, 2)
   e <- 2 * (sum(2 * y * pnorm(y) + 2 * dnorm(y)) - n / sqrt(pi) - mean(K * y))
   e
}
   
