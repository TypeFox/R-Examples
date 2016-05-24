###            Fit linear quantile models and linear quantile mixed models
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License or
#  any later version.
#
#  This program is distributed without any warranty,
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

".onAttach" <- function(lib, pkg) {
    if(interactive() || getOption("verbose"))
	packageStartupMessage(sprintf("Package %s (%s) loaded. Type citation(\"%s\") on how to cite this package\n", pkg,
		packageDescription(pkg)$Version, pkg))
}

##########################################################################################
# New generics

boot <- function(object, R = 50, seed = round(runif(1, 1, 10000)), startQR = FALSE) UseMethod("boot")
extractBoot <- function(object, which = "fixed") UseMethod("extractBoot")
predint <- function(object, level = 0, alpha = 0.05, R = 50, seed = round(runif(1, 1, 10000))) UseMethod("predint")

##########################################################################################
# Functions from contributed R packages (CRAN)

"is.positive.definite" <- function(m, tol, method = c("eigen", "chol")) 
{
# source package `corpcor' version 1.6.0 (Juliane Schaefer, Rainer Opgen-Rhein, Verena Zuber, A. Pedro Duarte Silva and Korbinian Strimmer)

    method = match.arg(method)
    if (!is.matrix(m)) 
        m = as.matrix(m)
    if (method == "eigen") {
        eval = eigen(m, only.values = TRUE)$values
        if (is.complex(eval)) {
            warning("Input matrix has complex eigenvalues!")
            return(FALSE)
        }
        if (missing(tol)) 
            tol = max(dim(m)) * max(abs(eval)) * .Machine$double.eps
        if (sum(eval > tol) == length(eval)) 
            return(TRUE)
        else return(FALSE)
    }
    if (method == "chol") {
        val = try(chol(m), silent = TRUE)
        if (class(val) == "try-error") 
            return(FALSE)
        else return(TRUE)
    }
}

"make.positive.definite" <- function(m, tol) 
{
# source package `corpcor' version 1.6.0 (Juliane Schaefer, Rainer Opgen-Rhein, Verena Zuber, A. Pedro Duarte Silva and Korbinian Strimmer)

    if (!is.matrix(m)) 
        m = as.matrix(m)
    d = dim(m)[1]
    if (dim(m)[2] != d) 
        stop("Input matrix is not square!")
    es = eigen(m)
    esv = es$values
    if (missing(tol)) 
        tol = d * max(abs(esv)) * .Machine$double.eps
    delta = 2 * tol
    tau = pmax(0, delta - esv)
    dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
    return(m + dm)
}

##

"permutations" <- function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) 
{

# source `gtools' version 2.6.2 (Gregory R. Warnes)

    if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
        0) 
        stop("bad value of n")
    if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
        0) 
        stop("bad value of r")
    if (!is.atomic(v) || length(v) < n) 
        stop("v is either non-atomic or too short")
    if ((r > n) & repeats.allowed == FALSE) 
        stop("r > n and repeats.allowed=FALSE")
    if (set) {
        v <- unique(sort(v))
        if (length(v) < n) 
            stop("too few different elements")
    }
    v0 <- vector(mode(v), 0)
    if (repeats.allowed) 
        sub <- function(n, r, v) {
            if (r == 1) 
                matrix(v, n, 1)
            else if (n == 1) 
                matrix(v, 1, r)
            else {
                inner <- Recall(n, r - 1, v)
                cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner), 
                  ncol = ncol(inner), nrow = nrow(inner) * n, 
                  byrow = TRUE))
            }
        }
    else sub <- function(n, r, v) {
        if (r == 1) 
            matrix(v, n, 1)
        else if (n == 1) 
            matrix(v, 1, r)
        else {
            X <- NULL
            for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n - 
                1, r - 1, v[-i])))
            X
        }
    }
    sub(n, r, v[1:n])
}

##

allVarsRec <- function(object)
{

# source `nlme' version 3.1-100 (Jose Pinheiro, Douglas Bates, Saikat DebRoy, Deepayan Sarkar and the R Development Core Team)

  if (is.list(object)) {
    unlist(lapply(object, allVarsRec))
  } else {
    all.vars(object)
  }
}

asOneFormula <- function(..., omit = c(".", "pi"))
{
# source `nlme' version 3.1-100 (Jose Pinheiro, Douglas Bates, Saikat DebRoy, Deepayan Sarkar and the R Development Core Team)

  names <- unique(allVarsRec((list(...))))
  names <- names[is.na(match(names, omit))]
  if (length(names)) {
    eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
  } else NULL
}

##

"gauss.quad" <- function(n, kind = "legendre", alpha = 0, beta = 0) 
{

# source `statmod' version 1.4.11 (Gordon Smyth)

    n <- as.integer(n)
    if (n < 0) 
        stop("need non-negative number of nodes")
    if (n == 0) 
        return(list(nodes = numeric(0), weights = numeric(0)))
    kind <- match.arg(kind, c("legendre", "chebyshev1", "chebyshev2", 
        "hermite", "jacobi", "laguerre"))
    i <- 1:n
    i1 <- i[-n]
    switch(kind, legendre = {
        muzero <- 2
        a <- rep(0, n)
        b <- i1/sqrt(4 * i1^2 - 1)
    }, chebyshev1 = {
        muzero <- pi
        a <- rep(0, n)
        b <- rep(0.5, n - 1)
        b[1] <- sqrt(0.5)
    }, chebyshev2 = {
        muzero <- pi/2
        a <- rep(0, n)
        b <- rep(0.5, n - 1)
    }, hermite = {
        muzero <- sqrt(pi)
        a <- rep(0, n)
        b <- sqrt(i1/2)
    }, jacobi = {
        ab <- alpha + beta
        muzero <- 2^(ab + 1) * gamma(alpha + 1) * gamma(beta + 
            1)/gamma(ab + 2)
        a <- i
        a[1] <- (beta - alpha)/(ab + 2)
        i2 <- 2:n
        abi <- ab + 2 * i2
        a[i2] <- (beta^2 - alpha^2)/(abi - 2)/abi
        b <- i1
        b[1] <- sqrt(4 * (alpha + 1) * (beta + 1)/(ab + 2)^2/(ab + 
            3))
        i2 <- i1[-1]
        abi <- ab + 2 * i2
        b[i2] <- sqrt(4 * i2 * (i2 + alpha) * (i2 + beta) * (i2 + 
            ab)/(abi^2 - 1)/abi^2)
    }, laguerre = {
        a <- 2 * i - 1 + alpha
        b <- sqrt(i1 * (i1 + alpha))
        muzero <- gamma(alpha + 1)
    })
    A <- rep(0, n * n)
    A[(n + 1) * (i - 1) + 1] <- a
    A[(n + 1) * (i1 - 1) + 2] <- b
    A[(n + 1) * i1] <- b
    dim(A) <- c(n, n)
    vd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(vd$vectors[1, ]))
    w <- muzero * w^2
    x <- rev(vd$values)
    list(nodes = x, weights = w)
}

"gauss.quad.prob" <- function(n, dist = "uniform", l = 0, u = 1, mu = 0, sigma = 1, 
    alpha = 1, beta = 1) 
{

# source `statmod' version 1.4.11 (Gordon Smyth)

    n <- as.integer(n)
    if (n < 0) 
        stop("need non-negative number of nodes")
    if (n == 0) 
        return(list(nodes = numeric(0), weights = numeric(0)))
    dist <- match.arg(dist, c("uniform", "beta1", "beta2", "normal", 
        "beta", "gamma"))
    if (n == 1) {
        switch(dist, uniform = {
            x <- (l + u)/2
        }, beta1 = , beta2 = , beta = {
            x <- alpha/(alpha + beta)
        }, normal = {
            x <- mu
        }, gamma = {
            x <- alpha * beta
        })
        return(list(nodes = x, weights = 1))
    }
    if (dist == "beta" && alpha == 0.5 && beta == 0.5) 
        dist <- "beta1"
    if (dist == "beta" && alpha == 1.5 && beta == 1.5) 
        dist <- "beta2"
    i <- 1:n
    i1 <- 1:(n - 1)
    switch(dist, uniform = {
        a <- rep(0, n)
        b <- i1/sqrt(4 * i1^2 - 1)
    }, beta1 = {
        a <- rep(0, n)
        b <- rep(0.5, n - 1)
        b[1] <- sqrt(0.5)
    }, beta2 = {
        a <- rep(0, n)
        b <- rep(0.5, n - 1)
    }, normal = {
        a <- rep(0, n)
        b <- sqrt(i1/2)
    }, beta = {
        ab <- alpha + beta
        a <- i
        a[1] <- (alpha - beta)/ab
        i2 <- 2:n
        abi <- ab - 2 + 2 * i2
        a[i2] <- ((alpha - 1)^2 - (beta - 1)^2)/(abi - 2)/abi
        b <- i1
        b[1] <- sqrt(4 * alpha * beta/ab^2/(ab + 1))
        i2 <- i1[-1]
        abi <- ab - 2 + 2 * i2
        b[i2] <- sqrt(4 * i2 * (i2 + alpha - 1) * (i2 + beta - 
            1) * (i2 + ab - 2)/(abi^2 - 1)/abi^2)
    }, gamma = {
        a <- 2 * i + alpha - 2
        b <- sqrt(i1 * (i1 + alpha - 1))
    })
    A <- rep(0, n * n)
    A[(n + 1) * (i - 1) + 1] <- a
    A[(n + 1) * (i1 - 1) + 2] <- b
    A[(n + 1) * i1] <- b
    dim(A) <- c(n, n)
    vd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(vd$vectors[1, ]))^2
    x <- rev(vd$values)
    switch(dist, uniform = x <- l + (u - l) * (x + 1)/2, beta1 = , 
        beta2 = , beta = x <- (x + 1)/2, normal = x <- mu + sqrt(2) * 
            sigma * x, gamma = x <- beta * x)
    list(nodes = x, weights = w)
}

"bandwidth.rq" <- function (p, n, hs = TRUE, alpha = 0.05) 
{

# source `quantreg' version 5.11 (Roger Koenker)

    x0 <- qnorm(p)
    f0 <- dnorm(x0)
    if (hs) 
        n^(-1/3) * qnorm(1 - alpha/2)^(2/3) * ((1.5 * f0^2)/(2 * 
            x0^2 + 1))^(1/3)
    else n^-0.2 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^0.2
}

##########################################################################################
# Asymmetric Laplace distribution

dal <- function(x, mu = 0, sigma = 1, tau = 0.5, log = FALSE) {

eps <- .Machine$double.eps^(2/3)
if(any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
if (tau == 0) tau <- eps
if (tau == 1) tau <- 1 - eps
if(sigma < 0) warning("Scale parameter 'sigma' is negative")

ind <- ifelse(x < mu, 1, 0)

val <- tau*(1-tau)/sigma * exp(-(x - mu)/sigma * (tau - ind))

if(log) log(val) else val

}

pal <- function(x, mu = 0, sigma = 1, tau = 0.5) {

eps <- .Machine$double.eps^(2/3)
if(any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
if (tau == 0) tau <- eps
if (tau == 1) tau <- 1 - eps
if(sigma < 0) warning("Scale parameter 'sigma' is negative")

ifelse(x < mu, tau * exp( (1 - tau) / sigma * (x - mu)),

1 - (1 - tau) *exp(- tau / sigma * (x - mu)))

}

ral <- function(n, mu = 0, sigma = 1, tau = 0.5) {

eps <- .Machine$double.eps^(2/3)
if(any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
if (tau == 0) tau <- eps
if (tau == 1) tau <- 1 - eps
if(sigma < 0) warning("Scale parameter 'sigma' is negative")

u <- runif(n)

x1 <- mu + sigma/(1-tau) * log(u/tau)

x2 <- mu - sigma/tau * log ((1-u) / (1-tau))

ifelse(x1 < mu, x1, x2)

}

qal <- function(x, mu = 0, sigma = 1, tau = 0.5) {

if (any(x > 1) | any(x < 0)) stop("x must be in [0,1]")

eps <- .Machine$double.eps^(2/3)
if(any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
if (tau == 0) tau <- eps
if (tau == 1) tau <- 1 - eps
if(sigma < 0) warning("Scale parameter 'sigma' is negative")

ifelse(x < tau, mu + (sigma/(1-tau))*log(x/tau),

mu - (sigma/tau)*log((1-x)/(1-tau)))

}

meanAL <- function(mu, sigma, tau){

eps <- .Machine$double.eps^(2/3)
if(any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
if (tau == 0) tau <- eps
if (tau == 1) tau <- 1 - eps
if(sigma < 0) warning("Scale parameter 'sigma' is negative")

mu + sigma*(1-2*tau)/(tau*(1-tau))

}

varAL <- function(sigma, tau) {

eps <- .Machine$double.eps^(2/3)
if(any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
if (any(tau == 0)) tau[tau == 0] <- eps
if (any(tau == 1)) tau[tau == 1] <- 1 - eps
if(sigma < 0) warning("Scale parameter 'sigma' is negative")

sigma^2*(1-2*tau+2*tau^2)/((1-tau)^2*tau^2)

}

invvarAL <- function(x, tau) {

eps <- .Machine$double.eps^(2/3)
if(any(tau > 1) | any(tau < 0)) stop("Parameter 'tau' must be in [0,1]")
if (tau == 0) tau <- eps
if (tau == 1) tau <- 1 - eps

sqrt(x*(tau*(1-tau))^2/(1-2*tau+2*tau^2))

}

mleAL <- function(x){

tau <- 0.5
m <- as.numeric(quantile(x, tau))
sigma <- 1

r <- 0
while(r < 1000){

m.last <- as.numeric(quantile(x, tau))
res <- x - m.last
sigma.last <- mean(res*(tau - ifelse(res<0,1,0)))
a <- mean(res*ifelse(res<=0,1,0))
tau.last <- (a + sqrt(a^2 - (mean(x) - m.last)*a))/(mean(x)-m.last)

dm <- abs(m-m.last)
ds <- abs(sigma-sigma.last)
dp <- abs(tau-tau.last)

if(all(c(dm,ds,dp)<0.0001)) break
	else {m <- m.last; tau <- tau.last; sigma <- sigma.last}
r <- r+1
}
list(m=m,sigma=sigma,tau=tau,r=r)
}
