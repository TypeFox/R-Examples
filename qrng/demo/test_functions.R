### Experimenting with test functions for MC ###################################

## Ideas generalized from Sobol, Asotsky (2003, "One more experiment on
## estimating high-dimensional integrals by quasi-Monte Carlo methods"),
## Mathematics and Computers in Simulation, 62 (3--6), 255--263

require(randtoolbox)
require(qrng)
require(copula)


### 1) Functions ###############################################################

### 1.1) QMC Marshall--Olkin sampling for Clayton and Gumbel ###################

##' @title Produces Samples of a d-dimensional Archimedean copula via
##'        Marshall--Olkin
##' @param u (n, d+1)-matrix sample of uniforms (first column mapped to V; we
##'        take the first one since this is typically best for QRNGs)
##' @param family Archimedean copula family
##' @param theta copula parameter
##' @return (n, d)-matrix of copula samples
##' @author Marius Hofert
rAC_MO <- function(u, family, theta)
{
    if(!is.matrix(u)) u <- rbind(u)
    dm <- dim(u)
    n <- dm[1]
    d <- dm[2] - 1
    V <- switch(family, # n-vector of frailties sampled via inversion
                "Clayton" = {
                    qgamma(u[,1], shape=1/theta)
                },
                "Gumbel" = {
                    ## takes long...
                    stabledist::qstable(u[,1], alpha=1/theta, beta=1, gamma=cos(pi/(2*theta))^theta,
                                        delta=theta==1, pm=1)
                },
                stop("'family' not implemented"))
    E <- qexp(u[,-1]) # (n, d)-matrix of Exp(1)
    (1 + E/matrix(rep(V, d), ncol=d))^(-1/theta) # (n, d)-matrix
}


### 1.2) Test functions ########################################################

## Note: They all have to integrate to 1

##' @title Test function (linear combination of powers; integrates to 1)
##' @param x (n, d) matrix of values in [0,1]
##' @param alpha d-vector of powers
##' @param w d-vector of weights, has to sum to 1
##' @return f(\bm{x}) = \sum_{j=1}^d w_j/(\alpha_j+1) x_j^{\alpha_j}
##' @author Marius Hofert
##' @note No checking because of performance
test_linear_power <- function(x, alpha, w) {
    n <- nrow(x)
    rowSums(rep(w*(alpha+1), each=n) * x^rep(alpha, each=n))
}

##' @title Test function (integrand K(C(U))+1/2; integrates to 1)
##' @param x (n, d) matrix of values in [0,1]
##' @param family Archimedean copula family
##' @param tau Kendall's tau
##' @return f(\bm{x}) = K(C(\bm{x})) + 1/2
##' @author Marius Hofert
##' @note No checking because of performance
test_Kendall <- function(x, family, tau) {
    theta <- iTau(getAcop(family), tau)
    d <- ncol(x)
    cop <- onacopulaL(family, nacList=list(theta, 1:d))
    pK(pCopula(x, copula=cop), cop=cop@copula, d=d) + 0.5
}

##' @title Test function (mapping copula data x to [0,1]^2 via Rosenblatt and
##'        then applying Faure's g_1)
##' @param x (n, d) matrix of copula values in [0,1]
##' @param alpha weights
##' @param family Archimedean copula family
##' @param tau Kendall's tau
##' @return f(\bm{u}) = \prod_{j=1}^d (|4u_j-2| + \alpha_j) / (1+\alpha_j) where
##'         u_j = C(x_j | x_1,\dots, x_{j-1}) (u = Rosenblatt-transformed x)
##' @author Marius Hofert
##' @note No checking because of performance
test_Faure <- function(x, alpha, family, tau) {
    d <- ncol(x)
    theta <- iTau(getAcop(family), tau)
    u <- rtrafo(x, cop=onacopulaL(family, nacList=list(theta, 1:d)), inverse=TRUE)
    alpha. <- matrix(rep(alpha, each=nrow(u)), ncol=d)
    apply((abs(4*u-2)+alpha.)/(1+alpha.), 1, prod)
}


### 1.3) Main function for testing the quality of (copula) quasi-random numbers

##' @title Absolute error when integrating a given test function
##' @param n *vector* of sample sizes
##' @param B number of bootstrap replications
##' @param d dimension
##' @param family Archimedean copula family
##' @param tau Kendall's tau (converted to the copula parameter)
##' @param test test function (a function in one argument which integrates to 1)
##' @param rng.method the random number generating method
##' @param sampling.method the copula sampling method (Marshall--Olkin or
##'        conditional distribution method)
##' @param ... additional arguments passed to sobol()
##' @return length(n)-vector of bootstrap-estimated mean absolute errors based
##'         on samples of size n
##' @author Marius Hofert
##' @note - Absolute error = relative error since test functions integrate to 1
##'       -
abs_err <- function(n, B, d, family=c("Clayton", "Gumbel"), tau, test,
                    rng.method=c("ghalton", "sobol", "prng"),
                    sampling.method=c("MO", "CDM"), ...)
{
    stopifnot(n >= 1, d >= 2, B >= 1, 0 < tau, tau < 1, is.function(test))
    family <- match.arg(family)
    rng.method <- match.arg(rng.method)
    sampling.method <- match.arg(sampling.method)

    ## Main
    res <- matrix(, nrow=B, ncol=length(n)) # (B, length(n))-matrix
    max.n <- max(n) # maximal n (for generating only one long sequence and 'picking out blocks')
    p <- if(sampling.method=="MO") d+1 else d
    for(b in 1:B) { # iterate over all bootstrap replications

        ## Generate (max(n), d)-matrix of U[0,1] samples
        U <- switch(rng.method,
                    "ghalton" = {
                        ghalton(max.n, d=p)
                    },
                    "sobol" = {
                        sobol(max.n, dim=p, scrambling=1)
                    },
                    "prng" = {
                        matrix(runif(max.n*p), ncol=p)
                    },
                    stop("Wrong 'rng.method'"))

        ## Build (max(n), d)-matrix of copula samples
        theta <- iTau(getAcop(family), tau) # convert tau to theta
        cop.data <- switch(sampling.method,
                           "MO" = {
                               rAC_MO(U, family=family, theta=theta)
                           },
                           "CDM" = {
                               cop <- onacopulaL(family, nacList=list(theta, 1:d))
                               rtrafo(U, cop=cop, inverse=TRUE)
                           },
                           stop("Wrong 'sampling.method'"))

        ## Apply the test function to each line of cop.data and compute the absolute error
        err <- test(cop.data) - 1 # max(n)-vector; 1 = true value

        ## Grab out mean absolute error (estimated based on different (sub-)sample sizes n)
        res[b,] <- vapply(seq_along(n), function(i) abs(mean(err[1:n[i]])),
                          numeric(1)) # length(n)-vector of estimated abs. error for bootstrap replication b

    }
    colMeans(res) # length(n)-vector of bootstrap-estimated mean absolute errors based on samples of size n
    ## Note: One could return two columns more here with bootstrapped confidence intervals...
}


### 2) Simulation and analysis #################################################

## Parameters
n <- seq(1e4, 2e5, by=1e4) # sample sizes
B <- 25 # (bootstrap) replications
rng <- c("ghalton", "sobol", "prng") # random number generators
sampling <- c("CDM", "MO") # sampling methods
family <- "Clayton" # Gumbel with MO takes too long due to numerical inversion of F for V
d <- c(2, 10) # dimensions
tau <- c(0.2, 0.5) # Kendall's tau
doPDF <- Sys.getenv("USER") == "mhofert" # logical indicating whether cropped PDF is produced

## Disclaimer
cat(paste0("Note: This simulation takes a couple of minutes.\n"))


### 2.1) Test function test_linear_power() #####################################

## Compute the absolute error (~ 11min)
filebasename <- paste0("sim_lin_pow_test_fun_",family)
file <- paste0(filebasename,".rda")
if(file.exists(file)) attach(file) else {
    err <- array(, dim=c(length(d), length(tau), length(rng), length(sampling), length(n)),
                 dimnames=list(d=d, tau=tau, rng=rng, sampling=sampling, n=n))
    cat(paste0("Starting simulation for test_linear_power():\n"))
    for(i in seq_along(d)) { # dimension d
        test_lin_pow <- function(x) test_linear_power(x, alpha=rep(2, d[i]), w=rep(1/d[i], d[i]))
        for(j in seq_along(tau)) { # Kendall's tau
            for(k in seq_along(rng)) { # rng.method
                for(l in seq_along(sampling)) { # sampling.method
                    set.seed(271)
                    rt <- system.time(err[i,j,k,l,] <-
                        abs_err(n, B=B, d=d[i], family=family, tau=tau[j],
                                test=test_lin_pow, rng.method=rng[k],
                                sampling.method=sampling[l]))
                    cat(paste0("d=",d[i],", tau=",tau[j],", rng=",rng[k],", sampling=",
                               sampling[l],": Done in ", round(rt[3]), "s\n"))
                }
            }
        }
    }
    save(err, file=file)
}

## Plot
for(i in seq_along(d)) { # dimension d
    for(j in seq_along(tau)) { # Kendall's tau
        mylty <- c("solid", "dashed", "dotted", "dotdash", "4C88C488", "12345678") # lty 1:4 + more distinguishable 2
        file <- paste0("fig_",filebasename,"_tau=",tau[j],"_d=",d[i],".pdf")
        if(doPDF) pdf(file=file, width=7, height=7)
        plot(n, err[i,j,1,1,], type="l", log="y", ylim=range(err[i,j,,,]),
             xlab=expression(italic(n)), ylab=substitute("Absolute error for"~italic(d)==D*","~tau==Tau,
                 list(D=d[i], Tau=tau[j])))
        lines(n, err[i,j,1,2,], lty=mylty[2])
        lines(n, err[i,j,2,1,], lty=mylty[3])
        lines(n, err[i,j,2,2,], lty=mylty[4])
        lines(n, err[i,j,3,1,], lty=mylty[5])
        lines(n, err[i,j,3,2,], lty=mylty[6])
        legend(130000, 6e-4, lty=mylty, bty="n",
               legend=c("G. Halton; CDM", "G. Halton; MO", "Sobol; CDM", "Sobol; MO",
               "Monte Carlo; CDM", "Monte Carlo; MO"))
        if(doPDF) simsalapar::dev.off.pdf(file=file)
    }
}


### 2.2) Test function test_Kendall() ##########################################

## Compute the absolute error
filebasename <- paste0("sim_Kendall_test_fun_",family)
file <- paste0(filebasename,".rda")
if(file.exists(file)) attach(file) else {
    err <- array(, dim=c(length(d), length(tau), length(rng), length(sampling), length(n)),
                 dimnames=list(d=d, tau=tau, rng=rng, sampling=sampling, n=n))
    cat(paste0("Starting simulation for test_Kendall():\n"))
    for(i in seq_along(d)) { # dimension d
        for(j in seq_along(tau)) { # Kendall's tau
            test_Ken <- function(x)
                test_Kendall(x, family=family, tau=tau[j])
            for(k in seq_along(rng)) { # rng.method
                for(l in seq_along(sampling)) { # sampling.method
                    set.seed(271)
                    rt <- system.time(err[i,j,k,l,] <-
                        abs_err(n, B=B, d=d[i], family=family, tau=tau[j],
                                test=test_Ken, rng.method=rng[k],
                                sampling.method=sampling[l]))
                    cat(paste0("d=",d[i],", tau=",tau[j],", rng=",rng[k],", sampling=",
                               sampling[l],": Done in ", round(rt[3]), "s\n"))
                }
            }
        }
    }
    save(err, file=file)
}

## Plot
for(i in seq_along(d)) { # dimension d
    for(j in seq_along(tau)) { # Kendall's tau
        mylty <- c("solid", "dashed", "dotted", "dotdash", "4C88C488", "12345678") # lty 1:4 + more distinguishable 2
        file <- paste0("fig_",filebasename,"_tau=",tau[j],"_d=",d[i],".pdf")
        if(doPDF) pdf(file=file, width=7, height=7)
        plot(n, err[i,j,1,1,], type="l", log="y", ylim=range(err[i,j,,,]),
             xlab=expression(italic(n)), ylab=substitute("Absolute error for"~italic(d)==D*","~tau==Tau,
                 list(D=d[i], Tau=tau[j])))
        lines(n, err[i,j,1,2,], lty=mylty[2])
        lines(n, err[i,j,2,1,], lty=mylty[3])
        lines(n, err[i,j,2,2,], lty=mylty[4])
        lines(n, err[i,j,3,1,], lty=mylty[5])
        lines(n, err[i,j,3,2,], lty=mylty[6])
        legend(130000, 3.5e-4, lty=mylty, bty="n",
               legend=c("G. Halton; CDM", "G. Halton; MO", "Sobol; CDM", "Sobol; MO",
               "Monte Carlo; CDM", "Monte Carlo; MO"))
        if(doPDF) simsalapar::dev.off.pdf(file=file)
    }
}


### 2.3) Test function test_Faure() ############################################

## Compute the absolute error
filebasename <- paste0("sim_Faure_test_fun_",family)
file <- paste0(filebasename,".rda")
if(file.exists(file)) attach(file) else {
    err <- array(, dim=c(length(d), length(tau), length(rng), length(sampling), length(n)),
                 dimnames=list(d=d, tau=tau, rng=rng, sampling=sampling, n=n))
    cat(paste0("Starting simulation for test_Faure():\n"))
    for(i in seq_along(d)) { # dimension d
        for(j in seq_along(tau)) { # Kendall's tau
            test_Fau <- function(x)
                test_Faure(x, alpha=rep(2, d[i]), family=family, tau=tau[j])
            for(k in seq_along(rng)) { # rng.method
                for(l in seq_along(sampling)) { # sampling.method
                    set.seed(271)
                    rt <- system.time(err[i,j,k,l,] <-
                        abs_err(n, B=B, d=d[i], family=family, tau=tau[j],
                                test=test_Ken, rng.method=rng[k],
                                sampling.method=sampling[l]))
                    cat(paste0("d=",d[i],", tau=",tau[j],", rng=",rng[k],", sampling=",
                               sampling[l],": Done in ", round(rt[3]), "s\n"))
                }
            }
        }
    }
    save(err, file=file)
}

## Plot
for(i in seq_along(d)) { # dimension d
    for(j in seq_along(tau)) { # Kendall's tau
        mylty <- c("solid", "dashed", "dotted", "dotdash", "4C88C488", "12345678") # lty 1:4 + more distinguishable 2
        file <- paste0("fig_",filebasename,"_tau=",tau[j],"_d=",d[i],".pdf")
        if(doPDF) pdf(file=file, width=7, height=7)
        plot(n, err[i,j,1,1,], type="l", log="y", ylim=range(err[i,j,,,]),
             xlab=expression(italic(n)), ylab=substitute("Absolute error for"~italic(d)==D*","~tau==Tau,
                 list(D=d[i], Tau=tau[j])))
        lines(n, err[i,j,1,2,], lty=mylty[2])
        lines(n, err[i,j,2,1,], lty=mylty[3])
        lines(n, err[i,j,2,2,], lty=mylty[4])
        lines(n, err[i,j,3,1,], lty=mylty[5])
        lines(n, err[i,j,3,2,], lty=mylty[6])
        legend(130000, 4e-4, lty=mylty, bty="n",
               legend=c("G. Halton; CDM", "G. Halton; MO", "Sobol; CDM", "Sobol; MO",
               "Monte Carlo; CDM", "Monte Carlo; MO"))
        if(doPDF) simsalapar::dev.off.pdf(file=file)
    }
}
