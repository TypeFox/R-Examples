# lr.R -- likelihood ratio tests for partially AR(1)
# Copyright (C) 2015 Matthew Clegg

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

loglik.par.fkf <- function (Y, rho, sigma_M, sigma_R, M0=0, R0=Y[1]) {
    # Given a sequence Y and a parameterization (rho, sigma_M, sigma_R) of an
    # associated PAR process, calculates the negative log likelihood that
    # Y would be observed under these process parameters.  

    if (length(Y) < 1) return(NA_real_)
    if (length(dim(Y)) > 0) Y <- Y[,1]
    Y <- coredata(Y)

    M0 <- as.numeric(M0)
    R0 <- as.numeric(R0)
    
    a0 <- c(M0, R0)
    dt <- matrix(0, 2, 1)
    ct <- matrix(0, 1, 1)
    Zt <- matrix(c(1, 1), 1, 2)
    GGt <- matrix(0, 1, 1)
    P0 <- matrix(c(sigma_M^2, 0, 0, sigma_R^2), 2, 2)
    Tt <- matrix(c(rho, 0, 0, 1), 2, 2)
    HHt <- matrix(c(sigma_M^2, 0, 0,sigma_R^2), 2, 2)

    sp <- list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, 
            GGt = GGt, HHt=HHt)
    
    ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
           Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = rbind(Y))
    
#    cat("Kalman gain = ", ans$Kt, "\n")         
    -ans$logLik 
}

loglik.par.ss <- function (Y, rho, sigma_M, sigma_R, M0=0, R0=Y[1]) {
    # Given a sequence Y and a parametrization (rho, sigma_M, sigma_R) of an 
    # associated PAR process, calculates the negative log likelihood that Y would
    # be observed under these process parameters, using a steady state 
    # Kalman filter.
    #
    # Despite the fact that this is hand-coded, the likelihood function that uses
    # the fast kalman filter implementation given above is about twice as fast,
    # even after turning on the byte compiler.  

    if (length(dim(Y)) > 0) Y <- Y[,1]
    Y <- coredata(Y)

    M0 <- as.numeric(M0)
    R0 <- as.numeric(R0)

    n <- length(Y)
    if (n < 1) return(NA_real_)
 
    K <- kalman.gain.par(rho, sigma_M, sigma_R)
    esumsq <- 0
    M <- M0
    R <- R0
    tvar <- sigma_M^2 + sigma_R^2
    
    for (i in 1:n) {
        xhat <- rho * M + R
        e <- Y[i] - xhat
        esumsq <- esumsq + e*e
        M <- rho * M + e * K[1]
        R <- R + e * K[2]
    }
    
    nll <- (n/2)*log(tvar * 2*pi) + esumsq/(2*tvar)
    
    nll
}

llst <- function(X, mu, sigma, df) {
    # Computes the negative log likelihood that (X - mu)/sigma is distributed
    # according to a t-distribution with df degrees of freedom.
    # Input values:
    #   X:     Vector of observations whose likelihood is to be computed
    #   mu:    The location parameter of the distribution
    #   sigma: The scale parameter -- similar to standard deviation
    #   df:    The degrees of freedom of the distribution
    #
    # This likelihood function has been optimized by computing the
    # values of the gamma function and some other constants only once,
    # rather than re-computing these values for each element of X.
    # This results in a factor of 10 (or more) speedup over the
    # naive implementation.
    
    if ((sigma <= 0) || (df <= 1) || (length(X) < 1)) return (NA_real_)

    # Following is the non-optimized version of the code:
    #  D <- dt((X - mu) / sigma, df) / sigma
    #  ll <- -sum(log(D))
    #  return(ll)

    # Following is an equivalent optimized version:
    N <- length(X)
    C <- N * (lgamma((df+1)/2) - 0.5*log(df * pi) - lgamma(df/2) - log(sigma))
    S <- -((df + 1)/2) * sum(log(1 + ((X-mu)/sigma)^2/df))
    ll <- -(S+C)
    ll
}

loglik.par.ss.t <- function (Y, rho, sigma_M, sigma_R, M0=0, R0=Y[1], nu=par.nu.default()) {
    # Given a sequence Y and a parametrization (rho, sigma_M, sigma_R) of an 
    # associated PAR process, calculates the negative log likelihood that Y would
    # be observed under these process parameters, using a steady state 
    # Kalman filter.
    
    if (length(dim(Y)) > 0) Y <- Y[,1]
    Y <- coredata(Y)

    M0 <- as.numeric(M0)
    R0 <- as.numeric(R0)

    K <- kalman.gain.par(rho, sigma_M, sigma_R)
    if (is.na(K[1])) return(NA_real_)
    
    n <- length(Y)
    if (n < 1) return(NA_real_)
    esumsq <- 0
    M <- M0
    R <- R0
    tvar <- sigma_M^2 + sigma_R^2
    tsd <- sqrt(tvar)
    
    evec <- numeric(n)
    for (i in 1:n) {
        xhat <- rho * M + R
        e <- Y[i] - xhat
        evec[i] <- e
        M <- rho * M + e * K[1]
        R <- R + e * K[2]
    }
    
#    nll <- llst(evec, 0, tsd * sqrt((nu - 2) / nu), nu)
    nll <- llst(evec, 0, tsd, nu)
    
    nll
}

loglik.par <- function (Y, rho, sigma_M, sigma_R, M0=0, R0=Y[1], 
    calc_method=c("css", "fkf", "ss", "sst", "csst"), nu=par.nu.default()) {
    # Given a sequence Y and a parameterization (rho, sigma_M, sigma_R) of an
    # associated PAR process, calculates the negative log likelihood that
    # Y would be observed under these process parameters.  The method used
    # for calculating the log likelihood is determined by "method":
    #   fkf:  Uses the Fast Kalman Filter (fkf) package
    #   ss:   Uses a steady state Kalman filter
    #   css:  Uses a steady state Kalman filter coded in C
    #   sst:  Uses a robust steady state Kalman filter, where the
    #         error terms are assumed to be t-distributed
    #   csst: Uses a robust steady state Kalman filter coded in C
    
    switch(match.arg(calc_method),
      fkf=loglik.par.fkf(Y, rho, sigma_M, sigma_R, M0, R0),
      ss=loglik.par.ss(Y, rho, sigma_M, sigma_R, M0, R0),
      css=loglik_par_c(Y, rho, sigma_M, sigma_R, M0, R0),
      sst=loglik.par.ss.t(Y, rho, sigma_M, sigma_R, M0, R0, nu),
      csst=loglik_par_t_c(Y, rho, sigma_M, sigma_R, M0, R0, nu))
    
}

likelihood_ratio.par <- function (X, # The series which is being fit
    robust=FALSE,            # If TRUE, robust estimations are performed                  
    null_model=c("rw", "ar1"),  # Specifies the null hypothesis
                             # rw = null model estimates sigma_R, assuming rho = sigma_M = 0.
                             #      This is the default.
                             # ar1 = null model estimates rho and sigma_M, assuming sigma_R = 0.
    opt_method=c("css", "fkf", "ss"),  
                             # Method to be used for calculating the log-likelihoods                         
                             #   css:  Uses a steady state Kalman filter coded in C
                             #   fkf:  Uses the Fast Kalman Filter (fkf) package
                             #   ss:   Uses a steady state Kalman filter
    nu=par.nu.default()      # If robust is TRUE, the degrees of freedom parameter                          
) {
    null_model <- match.arg(null_model)
    opt_method <- match.arg(opt_method)

    if (robust && opt_method == "fkf") stop("robust estimation not implemented for opt_method = fkf")
    f.alt <- fit.par(X, robust=robust, opt_method=opt_method, nu=nu)
    f.null <- fit.par(X, robust=robust, model=null_model, opt_method=opt_method, nu=nu)
    f.alt$negloglik - f.null$negloglik   
}

sample.likelihood_ratio.par <- function (n=500, rho=0.8, sigma_M=1.0, sigma_R=1.0, nrep=1000, 
    use.multicore=TRUE, robust=FALSE, nu=par.nu.default(), seed.start=0) {
    # Generates random sequences with the specified sample size, rho, sigma_M and sigma_R.
    # For each such random sequence, performs RW, AR(1) and PAR fits.  
    # Performs RW, AR(1) and PAR fits to random sequences with the specified sample size,
    # rho, sigma_M and sigma_R.  Creates a data.frame containing the statistics obtained
    # for each of the fits.

    param <- expand.grid(1:nrep, n, rho, sigma_M, sigma_R, robust, nu)
    colnames(param) <- c("rep", "n", "rho", "sigma_M", "sigma_R", "robust", "nu")
    param$seed <- (1:nrow(param)) + seed.start
    
    sample1 <- function (k) {
        set.seed(param$seed[k])
        X <- rpar(param$n[k], param$rho[k], param$sigma_M[k], param$sigma_R[k], 
            robust=param$robust[k], nu=param$nu[k])
        frw <- fit.par(X, robust=param$robust[k], model="rw", nu=param$nu[k])
        fmr <- fit.par(X, robust=param$robust[k], model="ar1", nu=param$nu[k])
        fpar <- fit.par(X, robust=param$robust[k], model="par", nu=param$nu[k])
        kpss <- suppressWarnings(kpss.test(X))
        c(param$n[k], param$rho[k], param$sigma_M[k], param$sigma_R[k], 
          param$robust[k], param$nu[k], param$seed[k],
          frw$rho, frw$sigma_M, frw$sigma_R, frw$negloglik,
          fmr$rho, fmr$sigma_M, fmr$sigma_R, fmr$negloglik,
          fpar$rho, fpar$sigma_M, fpar$sigma_R, fpar$negloglik, 
          fpar$negloglik - frw$negloglik, fpar$negloglik - fmr$negloglik,
          kpss$statistic, kpss$p.value, fpar$pvmr)
    }

    pvec_sample <- function(N) {
        lapply(N, function(k) sample1(k))
    }
    
    if (use.multicore) {
        samples <- parallel::pvec(1:nrow(param), pvec_sample)
    } else {
        samples <- lapply(1:nrow(param), sample1)
    }
    df <- as.data.frame(do.call("rbind", samples))
    colnames(df) <- c("n", "rho", "sigma_M", "sigma_R", "robust", "nu", "seed",
                      "rw_rho", "rw_sigma_M", "rw_sigma_R", "rw_negloglik",
                      "mr_rho", "mr_sigma_M", "mr_sigma_R", "mr_negloglik",  
                      "par_rho", "par_sigma_M", "par_sigma_R", "par_negloglik",
                      "rw_lrt", "mr_lrt", "kpss_stat", "kpss_p", "pvmr")
    df 
}

par.generate.likelihood_ratio.samples <- function (sample_dir = "samples",
    nrep = 10000,
    sample_size=c(50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 750, 800, 900, 1000, 1250, 
        1500, 1750, 2000, 2250, 2500), 
    rho=c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
    sigma_M=1.0,
    sigma_R=0.0,
    robust=c(FALSE, TRUE), 
    seed.start=1,
    skip.existing=FALSE,
    ...
) {
    # For each combination of sample size, rho and robust, generates
    # nrep likelihood_ratio samples with those parameters,
    # and writes the result to a file whose name is in the form
    #  sample_dir/PAR.SAMPLES.nn.rr[.ROB]
    params.df <- expand.grid(sample_size, robust, rho, sigma_M, sigma_R)
    colnames(params.df) <- c("sample_size", "robust", "rho", "sigma_M", "sigma_R")
    dir.create(sample_dir, showWarnings=FALSE, recursive=TRUE)
    for (i in 1:nrow(params.df)) {
        robstring <- if (params.df$robust[i]) ".ROB" else ""
        fn <- sprintf("%s/PAR.SAMPLES.%d.%02d.%02d.%02d%s", sample_dir, 
            params.df$sample_size[i], 
            floor(100.0 * params.df$rho[i]),
            floor(100.0 * params.df$sigma_M[i]),
            floor(100.0 * params.df$sigma_R[i]),
            robstring)
        if (!skip.existing || !file.exists(fn)) {
            printf("%s ", Sys.time())
            lrdf <- sample.likelihood_ratio.par(n=params.df$sample_size[i], rho=params.df$rho[i],
                sigma_M=params.df$sigma_M[i], sigma_R=params.df$sigma_R[i],
                nrep=nrep, robust=params.df$robust[i], seed.start=seed.start, ...)
            write.table(lrdf, fn)
            printf("%s\n", fn)
        }
        seed.start <- seed.start + nrep
    }
}

par.load.likelihood_ratio.samples <- function (sample_dir = "samples") {
    # Reads all of the samples in the specified samples dir, and returns
    # a data.frame containing the results
    
    files <- list.files(sample_dir, "PAR.SAMPLES", full.names=TRUE)
    SAMPLES <- do.call("rbind", lapply(files, function(f) read.table(f, header=TRUE)))
    SAMPLES$kpss_nstat <- -SAMPLES$kpss_stat
#    PAR.SAMPLES.DT <<- as.data.table(PAR.SAMPLES)
#    "Samples loaded into PAR.SAMPLES"
    if (exists("PAR.SAMPLES", envir=asNamespace("partialAR"))) {
        unlockBinding("PAR.SAMPLES", asNamespace("partialAR"))
        unlockBinding("PAR.SAMPLES.DT", asNamespace("partialAR"))
        rm(c("PAR.SAMPLES", "PAR.SAMPLES.DT"), envir=asNamespace("partialAR"))
    }
    G <- .GlobalEnv
    assign("PAR.SAMPLES", SAMPLES, envir=G)
    assign("PAR.SAMPLES.DT", data.table(SAMPLES), envir=G)
    invisible(SAMPLES)
}

par.joint.dist <- function (Lambda_R, Lambda_M, n=500, rho=0.9, robust=FALSE, 
    ar1test = c("lr", "kpss")) {
    # Given likelihood scores Lambda_R and Lambda_M, determines the proportion
    # of samples where the likelihood ratio for the random walk hypothesis is
    # less than Lambda_R and the likelihood ratio for the AR(1) hypothesis is
    # less than Lambda_M.  This function can be viewed as an inverse of 
    # par.find.joint.critical.values.
    rw_lrt <- mr_lrt <- kpss_nstat <- NULL
    
    ar1test <- match.arg(ar1test)
    pn <- n
    prho <- rho
    probust <- if (robust) 1 else 0

    if (!exists("PAR.SAMPLES.DT")) {
        if (!exists("PAR.SAMPLES")) PAR.SAMPLES <- par.load.likelihood_ratio.samples()
        PAR.SAMPLES.DT <- data.table(PAR.SAMPLES)
    }
    setkey(PAR.SAMPLES.DT, n, rho, robust)
    A <- PAR.SAMPLES.DT[list(pn, prho, probust),]
    if (nrow(A) == 0) stop ("No matching samples found")
    if (ar1test == "lr") {
        A[,sum(rw_lrt < Lambda_R & mr_lrt < Lambda_M)]/nrow(A)
    } else {
        A[,sum(rw_lrt < Lambda_R & kpss_nstat < Lambda_M)]/nrow(A)
    }        
#    A <- PAR.SAMPLES[PAR.SAMPLES$n == n & PAR.SAMPLES$rho==rho & PAR.SAMPLES$robust == robust,]
#    if (nrow(A) == 0) stop ("No matching samples found")
#    sum(A$rw_lrt < Lambda_R & A$mr_lrt < Lambda_M) / nrow(A)
}

par.find.joint.critical.values <- function (alpha, n=500, rho=0.9, robust=FALSE, nest=TRUE,
    ar1test = c("lr", "kpss")) {
    # Given an acceptance value alpha (or a vector of such acceptance values), 
    # calculates likelihood scores Lambda_R and Lambda_M such that the probability
    # is no more than alpha that a random sample will have likelihood ratio for the 
    # random walk no more than Lambda_R and likelihood ratio for the AR(1) hypothesis
    # no more than Lambda_M.  For the inverse, see par.joint.dist.
        
    ar1test <- match.arg(ar1test)
    pn <- n
    prho <- rho
    probust <- if (robust) 1 else 0
    rw <- NULL  # Make R CMD check happy
    mr <- NULL  # Make R CMD check happy

    if (!exists("PAR.SAMPLES.DT")) {
        if (!exists("PAR.SAMPLES")) PAR.SAMPLES <- par.load.likelihood_ratio.samples()
        PAR.SAMPLES.DT <- data.table(PAR.SAMPLES)
    }
    setkey(PAR.SAMPLES.DT, n, rho, robust)
    AM <- PAR.SAMPLES.DT[list(pn, prho, probust),]
    AR <- PAR.SAMPLES.DT[list(pn, 1, probust),]
    
    if (nrow(AM) == 0 || nrow(AR) == 0) stop ("No matching samples found") 

    if (ar1test == "lr") {
        AM <- AM[,c("rw_lrt", "mr_lrt"),with=FALSE]
        AR <- AR[,c("rw_lrt", "mr_lrt"),with=FALSE]
    } else {
        AM <- AM[,c("rw_lrt", "kpss_nstat"),with=FALSE]
        AR <- AR[,c("rw_lrt", "kpss_nstat"),with=FALSE]
    }
    setnames(AM, c("rw", "mr"))
    setnames(AR, c("rw", "mr"))
    
    LR_min <- -Inf
    LM_min <- -Inf
    
    find1 <- function(a) {
        cv <- function(p) {
            LR <- p[1]
            LM <- p[2]
            if (nest && ((LR < LR_min) || (LM < LM_min))) return(3)
            am <- AM[,sum(rw < LR & mr < LM)] / nrow(AM)
            ar <- AR[,sum(rw < LR & mr < LM)] / nrow(AR)
            (am - a)^2 + (ar - a)^2
        }
        
        LR0 <- max(quantile(AR$rw, a), LR_min)
        LM0 <- max(quantile(AM$mr, a), LM_min)
        p <- optim(c(LR0, LM0), cv)
        if (nest) {
            LR_min <<- p$par[1]
            LM_min <<- p$par[2]
        }
        QM <- AM[,sum(rw < p$par[1] & mr < p$par[2])] / nrow(AM)
        QR <- AR[,sum(rw < p$par[1] & mr < p$par[2])] / nrow(AR)
        v <- c(n, robust, a, LR0, LM0, p$par, QR, QM)
        names(v) <- c("n", "robust", "alpha", "C_R0", "C_M0", "C_R", "C_M", "QR", "QM")
        v
    }
    
    do.call("rbind", lapply(alpha, find1))
}

par.findall.joint.critical.values <- function (alpha=seq(0.01,0.99,by=0.01), rho=0.9, nest=TRUE,
    use.multicore=TRUE, ...) {
    if (!exists("PAR.SAMPLES")) PAR.SAMPLES <- par.load.likelihood_ratio.samples()
    idf <- expand.grid(n=sort(unique(PAR.SAMPLES$n)), robust=sort(unique(PAR.SAMPLES$robust)))
    if (use.multicore) lapply <- mclapply
    do.call("rbind", lapply(1:nrow(idf), function(k)
        par.find.joint.critical.values(alpha, idf$n[k], rho, idf$robust[k], nest=nest, ...)))
}

par.rw.pvalue <- function (
    stat.rw,                  # Statistic used for testing RW hypothesis.
                              # This is a likelihood ratio
    n,                        # Sample size
    robust=FALSE              # TRUE if robust models should be used
) {
    # Given a statistic on the RW hypothesis (e.g., the value of the
    # likelihood ratio function), returns the corresponding p-value.
    
    if (!robust) {
    	p.value <- quantile_table_interpolate(PAR.RWNULL.LRQT, n, stat.rw)
    } else {
    	p.value <- quantile_table_interpolate(PAR.RWNULL.ROB.LRQT, n, stat.rw)
    }

    p.value
}

par.mr.pvalue <- function (
    stat.mr,                  # Statistic used for testing AR(1) hypothesis
                              # If ar1test=="lr", then this is a likelihood ratio.
                              # If ar1test== "kpss", then this is the KPSS statistic
    n,                        # Sample size
    robust=FALSE,             # TRUE if robust models should be used
    ar1test=c("lr", "kpss")   # Test to be used against the null hypothesis of an AR(1) sequence
                              #  kpss = kpss.test, lr = likelihood ratio test
) {
    # Given a statistic on the AR(1) hypothesis, returns the corresponding 
    # p-value.  If ar1test == "lr", then the statistic stat.mr is the value
    # of the likelihood ratio.  If ar1test == "kpss", then the statistic stat.mr
    # is the value of the KPSS test.

    ar1test <- match.arg(ar1test)

    if (ar1test == "lr") {
        if (!robust) {
        	p.value <- quantile_table_interpolate(PAR.MRNULL.LRQT, n, stat.mr)
        } else {
        	p.value <- quantile_table_interpolate(PAR.MRNULL.ROB.LRQT, n, stat.mr)
        }
    } else {
        if (!robust) {
        	p.value <- quantile_table_interpolate(PAR.MRNULL.KPSS, n, stat.mr)
        } else {
        	p.value <- quantile_table_interpolate(PAR.MRNULL.ROB.KPSS, n, stat.mr)
        }
    }
    
    p.value
}

par.joint.pvalue <- function (
    stat.rw,                  # Statistic used for testing RW hypothesis.
                              # This is a likelihood ratio
    stat.mr,                  # Statistic used for testing AR(1) hypothesis
                              # If ar1test=="lr", then this is a likelihood ratio.
                              # If ar1test== "kpss", then this is the KPSS statistic
    n,                        # Sample size
    robust=FALSE,             # TRUE if robust models should be used
    ar1test=c("lr", "kpss")   # Test to be used against the null hypothesis of an AR(1) sequence
                              #  kpss = kpss.test, lr = likelihood ratio test
) {
    # Given statistics on the RW hypothesis and AR(1) hypothesis, calculates
    # a p-value associated with those statistics.  The p-value is an estimate
    # of the probability that a random sequence satisfying the null hypothesis
    # of either RW or AR(1) would have statistics more extreme than the
    # statistics (stat.rw, stat.mr)

    ar1test <- match.arg(ar1test)

    if (!exists("PAR.JOINT.CRITICAL.VALUES")) stop ("Could not find table of critical values")
    
    if (ar1test == "lr") {
        AJCV <- PAR.JOINT.CRITICAL.VALUES
    } else {
        AJCV <- PAR.JOINT.CRITICAL.VALUES.KPSS
    }

    nvals <- unique(AJCV[,"n"])
    probust <- if (robust) 1 else 0
    if (!any(nvals <= n)) {
        warning("Sample size too small (", n, ") to provide accurate p-value")
        nlower <- 1
        plower <- 1
    } else {
        nlower <- max(nvals[nvals <= n])
        lower.matches <- (AJCV[,"n"] == nlower) & 
            (AJCV[,"robust"] == probust) &
            (stat.rw <= AJCV[,"C_R"]) &
            (stat.mr <= AJCV[,"C_M"])
        if (!any(lower.matches)) {
            plower <- 1
        } else {
            plower <- min(AJCV[lower.matches,"alpha"])
        }
    }
    
    if (any(nvals >= n)) {
        nupper <- min(nvals[nvals > n])
        upper.matches <- (AJCV[,"n"] == nlower) & 
            (AJCV[,"robust"] == probust) &
            (stat.rw <= AJCV[,"C_R"]) &
            (stat.mr <= AJCV[,"C_M"])
        if (!any(upper.matches)) {
            pupper <- 1
        } else {
            pupper <- min(AJCV[upper.matches,"alpha"])
        }
    } else {
        nupper <- n
        pupper <- plower
    }
    
    pval.joint <- (nupper - n)/(nupper - nlower) * plower +
        (n - nlower) / (nupper - nlower) * pupper

    pval.joint
}

test.par.nullrw <- function (Y, robust=FALSE) {
    # Tests the hypothesis that Y is a random walk vs. the alternative
    # that Y is an PAR series.
    
    DNAME <- deparse(substitute(Y))
    STAT <- likelihood_ratio.par (Y, robust=robust, null_model="rw")
    names(STAT) <- "LL"
    PVAL <- par.rw.pvalue(STAT, length(Y), robust)
    METHOD <- ifelse(robust,
        "Likelihood ratio test of Robust Random Walk vs Robust Almost AR(1)",
        "Likelihood ratio test of Random Walk vs Almost AR(1)")
    alternative <- ifelse(robust, "Robust Almost AR(1)", "Almost AR(1)")

    structure(list(statistic = STAT, alternative = alternative, 
        p.value = PVAL, method = METHOD, data.name = DNAME),
        class = "htest")
}

test.par.nullmr <- function (Y, robust=FALSE, ar1test=c("lr", "kpss")) {
    # Tests whether a sequence Y conforms to the PAR model.  The null
    # hypothesis is that Y is AR(1).  If ar1test = "kpss", then the KPSS
    # test is used.  Otherwise, the likelihood ratio test is used.

    ar1test <- match.arg(ar1test)
    rstat <- paste(ar1test, ifelse(robust, "r", ""), sep="")
    
    DNAME <- deparse(substitute(Y))
    STATVAL <- switch(ar1test,
        lr = likelihood_ratio.par(Y, robust=robust, null_model="ar1"),
        kpss = suppressWarnings(-kpss.test(Y)$statistic))
    names(STATVAL) <- switch(ar1test, lr = "LL", kpss = "KPSS")
    METHOD <- switch (rstat,
        lr = "Likelihood ratio test of AR(1) vs Almost AR(1)",
        lrr = "Likelihood ratio test of Robust AR(1) vs Robust Almost AR(1)",
        kpss = "KPSS test of AR(1) vs Almost AR(1)",
        kpssr = "KPSS test of Robust AR(1) vs Robust Almost AR(1)")
    alternative <- ifelse (robust, "Almost AR(1)", "Robust Almost AR(1)")
    PVAL <- par.mr.pvalue(STATVAL, length(Y), robust=robust, ar1test=ar1test)
        
    structure(list(statistic = STATVAL, alternative = alternative, 
        p.value = PVAL, method = METHOD, data.name = DNAME),
        class = "htest")
}

test.par <- function (Y, 
    alpha=0.05,               # Critical value to be used in deciding whether to reject the null.
    null_hyp=c("rw", "ar1"),  # Specifies the null hypothesis.  Can be either or both.
                              # rw = null model estimates sigma_R, assuming rho = sigma_M = 0.
                              # ar1 = null model estimates rho and sigma_M, assuming sigma_R = 0.
                              # Default is both
    ar1test=c("lr", "kpss"),  # Test to be used against the null hypothesis of an AR(1) sequence
                              #  kpss = kpss.test, lr = likelihood ratio test
    robust=FALSE              # TRUE if robust models should be used
) {
    # Tests the null hypothesis that Y is either a random walk or an AR(1)
    # series (or both) against the alternative hypothesis that Y is PAR.

    DNAME <- deparse(substitute(Y))
    ar1test <- match.arg(ar1test)
    
    if (is(Y,"par.fit")) {
        Yfit <- Y
        Y <- Yfit$data
        if (missing(robust)) robust <- Yfit$robust
    }

    if (length(Y) < 50) stop("Input data is of insufficient length")
    
    if (length(null_hyp) == 1) {
        if (null_hyp == "rw") return(test.par.nullrw(Y, robust=robust))
        if ((null_hyp == "ar1") || (null_hyp == "mr")) return(test.par.nullmr(Y, ar1test=ar1test, robust=robust))
        stop("Error in parameter null_hyp: ", null_hyp)
    }
    
    frw <- fit.par(Y, model="rw", robust=robust)
    fpar <- fit.par(Y, model="par", rho.max=1, robust=robust)
    stat.rw <- fpar$negloglik - frw$negloglik
        
    if (ar1test == "lr") {
        fmr <- fit.par(Y, model="ar1", rho.max=1, robust=robust)
        stat.mr <- fpar$negloglik - fmr$negloglik
    } else {
        stat.mr <- suppressWarnings(-kpss.test(Y)$statistic)
    }
    STAT <- c(stat.rw, stat.mr)
    
    pval.joint <- par.joint.pvalue(stat.rw, stat.mr, length(Y), 
        ar1test=ar1test, robust=robust)
        
    if (!robust) {
    	pval.rw <- quantile_table_interpolate(PAR.RWNULL.LRQT, length(Y), stat.rw)
        if (ar1test == "lr") {
            pval.mr <- quantile_table_interpolate(PAR.MRNULL.LRQT, length(Y), stat.mr)
        } else {
            pval.mr <- quantile_table_interpolate(PAR.MRNULL.KPSS, length(Y), stat.mr)
        }
        METHOD <- "Test of [Random Walk or AR(1)] vs Almost AR(1)"
        alternative <- "Almost AR(1)"
    } else {
    	pval.rw <- quantile_table_interpolate(PAR.RWNULL.ROB.LRQT, length(Y), stat.rw)
        if (ar1test == "lr") {
            pval.mr <- quantile_table_interpolate(PAR.MRNULL.ROB.LRQT, length(Y), stat.mr)
        } else {
            pval.mr <- quantile_table_interpolate(PAR.MRNULL.ROB.KPSS, length(Y), stat.mr)
        }
        METHOD <- "Test of [Robust Random Walk or Robust AR(1)] vs Robust Almost AR(1)"
        alternative <- "Robust Almost AR(1)"
    }
    
    if (ar1test == "lr") {
        METHOD <- paste(METHOD, "[LR test for AR1]")
    } else {
        METHOD <- paste(METHOD, "[KPSS test for AR1]")
    }
    
    PVAL <- c(RW=pval.rw, AR1=pval.mr, PAR=pval.joint)
    
    names(STAT) <- "LL"
    structure(list(statistic = STAT, alternative = alternative, 
        p.value = PVAL, method = METHOD, data.name = DNAME, alpha=alpha, robust=robust),
        class = c("partest", "htest"))
}

print.partest <- function (x, ...) {
    # See stats:::print.htest
    print.partest.internal (x, ...)
}

print.partest.internal <- function (AT, alpha) {
    # See stats:::print.htest
    cat("\n")
    cat(strwrap(AT$method, prefix = "\t"), sep = "\n")
    cat("\n")
    cat("data:  ", AT$data.name, "\n", sep = "")
    cat("\n")
    
    if (AT$robust) {
        h0a <- "Robust RW"
        h0b <- "Robust AR(1)"
        h1 <- "Robust PAR"
    } else {
        h0a <- "Random Walk"
        h0b <- "AR(1)"
        h1 <- "PAR"
    }
    
    cat(sprintf("%-12s %20s %10s\n", "Hypothesis", "Statistic", "p-value"))
    cat(sprintf("%-12s %20.2f %10.3f\n", h0a, AT$statistic[1], AT$p.value[1]))
    cat(sprintf("%-12s %20.2f %10.3f\n", h0b, AT$statistic[2], AT$p.value[2]))
    cat(sprintf("%-12s %20s %10.3f\n", "Combined", "", AT$p.value[3]))
    
#    cat("\n")
#    cat("Acceptance Regions:\n")
#    if (AT$p.value[1] < AT$p.value[2]) {
#        cat(sprintf("  Accept hypothesis of %12s: alpha < %.3f\n", h0a, AT$p.value[1]))
#        cat(sprintf("  Accept hypothesis of %12s:         %.3f < alpha < %.3f\n", 
#            h0b, AT$p.value[1], AT$p.value[2]))
#        cat(sprintf("  Accept hypothesis of %12s:                         %.3f < alpha\n", h1, AT$p.value[2]))
#    } else {
#        cat(sprintf("  Accept hypothesis of %12s: alpha < %.3f\n", h0b, AT$p.value[2]))
#        cat(sprintf("  Accept hypothesis of %12s:         %.3f < alpha < %.3f\n", 
#            h0a, AT$p.value[2], AT$p.value[1]))
#        cat(sprintf("  Accept hypothesis of %12s:                         %.3f < alpha\n", h1, AT$p.value[1]))
#    }

#    if (!is.na(AT$alpha)) {
#        cat("\n")
#        cat(sprintf("For alpha = %.3f, the accepted hypothesis is: \n  ", AT$alpha))
#        cat(which.hypothesis.partest (AT))
#    }
    
    cat("\n")
    invisible(AT)
}

which.hypothesis.partest <- function(AT) {
    if (!is(AT, "partest")) {
        stop("which.hypothesis.partest must be called with an object of class partest")
    }
    
    if (is.na(AT$alpha)) {
        return(NA)
    } 

    if (AT$p.value[3] < AT$alpha) {
        result <- "PAR"
#    } else if (AT$p.value[1] < AT$p.value[2]) {
    } else if (AT$p.value[1] > AT$alpha) {
        result <- "RW"
    } else if (AT$p.value[2] > AT$alpha) {
        result <- "AR1"
    } else if (AT$p.value[1] > AT$p.value[2]) {
        result <- "RW"
    } else {
        result <- "AR1"
    }
    
#    if (max(AT$p.value[1], AT$p.value[2]) < AT$alpha) {
#        result <- "PAR"
#    } else if (AT$alpha < AT$p.value[1]) {
#        result <- "RW"
#    } else {
#        result <- "AR1"
#    }
    
    if (AT$robust) result <- paste("R", result, sep="")
    
    result
}

print.par.lrt <- function (latex = FALSE, robust=FALSE) {
      # Pretty prints a table of critical values for the likelihood
      # ratio tests.  If latex is TRUE, then also outputs latex formatting
      # commands.
      

    crit_values<-c(0.01, 0.05, 0.10)
    sample_sizes<-c(50, 100, 250, 500, 1000, 2500)
              
    ctext <- function (str, n) {
        if (nchar(str) == n) return(str)
        if (nchar(str) > n) return(str, 1, n)
        
        nleft <- floor((n - nchar(str))/2)
        nright <- n - (nchar(str) + nleft)
        left_pad <- paste(rep(" ", nleft), collapse="")
        right_pad <- paste(rep(" ", nright), collapse="")
        pstr <- paste(left_pad, str, right_pad, sep="")
        pstr
    }
    
    if (latex) {
        cat("\\begin{table}\n")
        cat("\\begin{tabular}{crrr|rrr}\n")
        cat(" & \\multicolumn{3}{c}{NULL: Random Walk} & \\multicolumn{3}{c}{NULL: AR(1)} \\\\ \n")
    } else {
        cat(ctext("Critical Values for Likelihood Ratio Tests", 80), "\n")
        cat(ctext("Single Hypothesis Test", 80), "\n")
        if (robust) cat(ctext("Robust Model", 80), "\n")
        cat("\n")
        cat(sprintf("%10s %30s | %30s\n", "", ctext("NULL: Random Walk", 30), 
            ctext("NULL: AR(1)", 30)))
    }
    
    if (latex) {
        cat(" ")
        for (p in crit_values) {
            cat(sprintf(" & \\multicolumn{1}{c}{p=%.2f}", p))
        }
        for (p in crit_values) {
            cat(sprintf(" & p=%.2f", p))
        }
        cat("\\\\ \n")
    } else {
        cat(sprintf("%10s ", ""))
        for (p in crit_values) {
            pstr <- sprintf("p=%.2f", p)
            cat(sprintf("%9s ", pstr))
        }
        cat(" | ")
        for (p in crit_values) {
            pstr <- sprintf("p=%.2f", p)
            cat(sprintf("%9s ", pstr))
        }
        cat("\n")
    }

    if (latex) {
        cat("\\hline \n")
    } else {
        cat(sprintf ("%10s   %s\n", "", paste(rep("-", 60), collapse="")))
    }

    if (robust) {
        MR0 <- PAR.RWNULL.ROB.LRQT
        RW0 <- PAR.MRNULL.ROB.LRQT
    } else {
        MR0 <- PAR.RWNULL.LRQT
        RW0 <- PAR.MRNULL.LRQT
    }
    
    for (n in sample_sizes) {
        nstr <- sprintf("n=%d", n)
        if (latex) {
            cat(nstr, " ")
        } else {
            cat(sprintf("%10s ", nstr))
        }
        j <- which(MR0[1,] == n)
        for (p in crit_values) {
            i <- which(MR0[,1] == p)
            if (latex) cat (" & ")
            cat(sprintf("%9.1f ", MR0[i,j]))
        }
        if (!latex) cat(" | ")
        j <- which(RW0[1,] == n)
        for (p in crit_values) {
            i <- which(RW0[,1] == p)
            if (latex) cat (" & ")
            cat(sprintf("%9.1f ", RW0[i,j]))
        }
        if (latex) cat (" \\\\ ")
        cat("\n")
    }  
    
    if (latex) {
        cat("\\end{tabular}\n") 
        cat("\\caption{Critical Values for Likelihood Ratio Tests} \n")
        cat("\\caption*{For each sample size, 40,000 random walks were generated, and then the\n")
        cat("likelihood ratios were calculated under the hypothesis of a random walk\n")
        cat("(left panel) and under the hypothesis of an AR(1) series (right panel).\n")
        cat("For the hypothesis of an AR(1) series, it was found that the critical values\n")
        cat("depend upon the value of $\\rho$, and that as $\\rho$ increases, the critical values\n")
        cat("for a given quantile decrease.  Thus, by using the limiting case of a random walk\n")
        cat("when computing critical values for the AR(1) case, a conservative estimate is\n")
        cat("obtained.}\n")
        cat("\\end{table}")
    }        
}

print.par.lrt.rw <- function (latex = FALSE, robust=FALSE) {
      # Pretty prints a table of critical values for the likelihood
      # ratio test when the null hypothesis is a random walk.  
      # If latex is TRUE, then also outputs latex formatting
      # commands.
      
    crit_values<-c(0.01, 0.05, 0.10)
    sample_sizes<-c(50, 100, 250, 500, 1000, 2500)
              
    ctext <- function (str, n) {
        if (nchar(str) == n) return(str)
        if (nchar(str) > n) return(str, 1, n)
        
        nleft <- floor((n - nchar(str))/2)
        nright <- n - (nchar(str) + nleft)
        left_pad <- paste(rep(" ", nleft), collapse="")
        right_pad <- paste(rep(" ", nright), collapse="")
        pstr <- paste(left_pad, str, right_pad, sep="")
        pstr
    }
    
    if (latex) {
        cat("\\begin{tabular}{crrr}\n")
        cat(" & \\multicolumn{3}{c}{NULL: Random Walk} \\\\ \n")
    } else {
        cat(ctext("Critical Values for Likelihood Ratio Tests", 50), "\n")
        if (robust) cat(ctext("Robust Model", 50), "\n")
        cat(ctext("Null hypothesis:  Random Walk", 50), "\n\n")
    }
    
    if (latex) {
        cat(" ")
        for (p in crit_values) {
            cat(sprintf(" & \\multicolumn{1}{c}{p=%.2f}", p))
        }
        cat("\\\\ \n")
    } else {
        cat(sprintf("%10s ", ""))
        for (p in crit_values) {
            pstr <- sprintf("p=%.2f", p)
            cat(sprintf("%9s ", pstr))
        }
        cat("\n")
    }

    if (latex) {
        cat("\\hline \n")
    } else {
        cat(sprintf ("%10s   %s\n", "", paste(rep("-", 28), collapse="")))
    }

    if (robust) {
        MR0 <- PAR.RWNULL.ROB.LRQT
    } else {
        MR0 <- PAR.RWNULL.LRQT
    }
     
    for (n in sample_sizes) {
        nstr <- sprintf("n=%d", n)
        if (latex) {
            cat(nstr, " ")
        } else {
            cat(sprintf("%10s ", nstr))
        }
        j <- which(MR0[1,] == n)
        for (p in crit_values) {
            i <- which(MR0[,1] == p)
            if (latex) cat (" & ")
            cat(sprintf("%9.1f ", MR0[i,j]))
        }
        if (latex) cat (" \\\\ ")
        cat("\n")
    }  
    
    if (latex) {
        cat("\\end{tabular}\n") 
    }        
}

print.par.lrt.mr <- function (latex = FALSE, robust=FALSE) {
      # Pretty prints a table of critical values for the likelihood
      # ratio test when the null hypothesis is AR(1).  
      # If latex is TRUE, then also outputs latex formatting
      # commands.
      
    crit_values<-c(0.01, 0.025, 0.05)
    sample_sizes<-c(50, 100, 250, 500, 1000, 2500)
    rho_values <- c(0.50, 0.80, 0.90)
              
    
    if (latex) {
        cat("\\begin{tabular}{crrr|rrr|rrr}\n")
        for (r in rho_values) {
            printf(" & \\multicolumn{3}{c}{\\rho = %.2f} ", r)
        } 
        cat(" \\\\ \n")
    } else {
        cat(ctext("Critical Values for Likelihood Ratio Tests", 120), "\n")
        cat(ctext("Null Hypothesis: AR(1)", 120), "\n")
        if (robust) cat(ctext("Robust Model", 120), "\n")
        cat("\n")
        printf("%10s ", "")
        for (r in rho_values) {
            printf("  %30s ", ctext(sprintf("rho = %.2f", r), 30))
        }
        printf("\n")
    }
    
    if (latex) {
        cat(" ")
        for (r in rho_values) {
            for (p in crit_values) {
                printf(" & \\multicolumn{1}{c}{p=%.2f}", p)
            }
        }
        cat("\\\\ \n")
    } else {
        cat(sprintf("%10s ", ""))
        for (r in rho_values) {
            for (p in crit_values) {
                pstr <- sprintf("p=%.3f", p)
                printf("%9s ", pstr)
            }
            printf("   ")
        }
        cat("\n")
    }

    if (latex) {
        cat("\\hline \n")
    } else {
        printf ("%10s   %s\n", "", paste(rep("-", 95), collapse=""))
    }

    print_values <- function (A, size) {
        j <- which(A[1,] == size)
        for (p in crit_values) {
            i <- which(A[,1] == p)
            if (latex) cat (" & ")
            cat(sprintf("%9.2f ", A[i,j]))
        }
    }

    if (!exists("PAR.SAMPLES")) PAR.SAMPLES <- par.load.likelihood_ratio.samples()
    
    r <- if (robust) 1 else 0

    MR50 <- quantile.table.from.samples("mr_lrt", 
        PAR.SAMPLES[PAR.SAMPLES$rho==0.5 & PAR.SAMPLES$robust == r,],
        quantiles = c(0.01, 0.025, 0.05))
    MR80 <- quantile.table.from.samples("mr_lrt", 
        PAR.SAMPLES[PAR.SAMPLES$rho==0.8 & PAR.SAMPLES$robust == r,],
        quantiles = c(0.01, 0.025, 0.05))
    MR90 <- quantile.table.from.samples("mr_lrt", 
        PAR.SAMPLES[PAR.SAMPLES$rho==0.9 & PAR.SAMPLES$robust == r,],
        quantiles = c(0.01, 0.025, 0.05))
    
    for (n in sample_sizes) {
        nstr <- sprintf("n=%d", n)
        if (latex) {
            cat(nstr, " ")
        } else {
            printf("%10s ", nstr)
        }
        print_values(MR50, n)
        if (!latex) cat(" | ")
        print_values(MR80, n)
        if (!latex) cat(" | ")
        print_values(MR90, n)
        
        cat("\n")
    }  
    
    if (latex) {
        cat("\\end{tabular}\n") 
    }        
}

print.par.kpss.mr <- function (latex = FALSE, robust=FALSE) {
      # Pretty prints a table of critical values for the KPSS
      # unit root test when the null hypothesis is AR(1).  
      # If latex is TRUE, then also outputs latex formatting
      # commands.
      
#    crit_values<-c(0.01, 0.025, 0.05)
    crit_values<-c(0.01, 0.05, 0.1)
    sample_sizes<-c(50, 100, 250, 500, 1000, 2500)
    rho_values <- c(0.50, 0.80, 0.90)            
    
    if (latex) {
        cat("\\begin{tabular}{crrr|rrr|rrr}\n")
        for (r in rho_values) {
            printf(" & \\multicolumn{3}{c}{\\rho = %.2f} ", r)
        } 
        cat(" \\\\ \n")
    } else {
        cat(ctext("Critical Values for KPSS Test", 120), "\n")
        cat(ctext("Null Hypothesis: AR(1)", 120), "\n")
        if (robust) cat(ctext("Robust Model", 120), "\n")
        cat("\n")
        printf("%10s ", "")
        for (r in rho_values) {
            printf("  %30s ", ctext(sprintf("rho = %.2f", r), 30))
        }
        printf("\n")
    }
    
    if (latex) {
        cat(" ")
        for (r in rho_values) {
            for (p in crit_values) {
                printf(" & \\multicolumn{1}{c}{p=%.2f}", p)
            }
        }
        cat("\\\\ \n")
    } else {
        cat(sprintf("%10s ", ""))
        for (r in rho_values) {
            for (p in crit_values) {
                pstr <- sprintf("p=%.3f", p)
                printf("%9s ", pstr)
            }
            printf("   ")
        }
        cat("\n")
    }

    if (latex) {
        cat("\\hline \n")
    } else {
        printf ("%10s   %s\n", "", paste(rep("-", 95), collapse=""))
    }

    print_values <- function (A, size) {
        j <- which(A[1,] == size)
        for (p in crit_values) {
            i <- which(A[,1] == p)
            if (latex) cat (" & ")
            cat(sprintf("%9.2f ", A[i,j]))
        }
    }

    if (!exists("PAR.SAMPLES")) PAR.SAMPLES <- par.load.likelihood_ratio.samples()
    
    r <- if (robust) 1 else 0

    MR50 <- quantile.table.from.samples("kpss_nstat", 
        PAR.SAMPLES[PAR.SAMPLES$rho==0.5 & PAR.SAMPLES$robust == r,],
        quantiles = crit_values)
    MR80 <- quantile.table.from.samples("kpss_nstat", 
        PAR.SAMPLES[PAR.SAMPLES$rho==0.8 & PAR.SAMPLES$robust == r,],
        quantiles = crit_values)
    MR90 <- quantile.table.from.samples("kpss_nstat", 
        PAR.SAMPLES[PAR.SAMPLES$rho==0.9 & PAR.SAMPLES$robust == r,],
        quantiles = crit_values)
    
    for (n in sample_sizes) {
        nstr <- sprintf("n=%d", n)
        if (latex) {
            cat(nstr, " ")
        } else {
            printf("%10s ", nstr)
        }
        print_values(MR50, n)
        if (!latex) cat(" | ")
        print_values(MR80, n)
        if (!latex) cat(" | ")
        print_values(MR90, n)
        
        cat("\n")
    }  
    
    if (latex) {
        cat("\\end{tabular}\n") 
    }        
}

print.par.lrt.quantile_quadrants <- function (latex = FALSE, robust=FALSE, rho = 0.9) {
    # Constructs a 2 x 2 table showing the quantiles of the likelihood
    # ratio function when the samples are either (a) generated from a
    # random walk process, or (b) generated from an AR(1) process with
    # mean-reversion parameter rho.  For each of the two groups of samples,
    # computes the likelihood ratio for a random walk null and an AR(1) null.
    # Then tabulates the quantiles of the likelihood ratio scores found in
    # this manner.

    q <- c(0.01, 0.05, 0.10)
    rob <- if (robust) 1 else 0
    sample_sizes<-c(50, 100, 250, 500, 1000, 2500)
    if (!exists("PAR.SAMPLES")) PAR.SAMPLES <- par.load.likelihood_ratio.samples()

    print_values <- function (A, size) {
        j <- which(A[1,] == size)
        for (p in q) {
            i <- which(A[,1] == p)
            if (latex) cat (" & ")
            printf("%9.2f ", A[i,j])
        }
    }
    
    print.two.panels <- function (r) {
        # Prints two panels of the likelihood ratio quadrants, corresponding
        # to the results obtained when the samples are generated from a
        # process where rho = r (r = 1 for random walk).
        
        RW <- quantile.table.from.samples("rw_lrt",
            PAR.SAMPLES[PAR.SAMPLES$rho==r & PAR.SAMPLES$robust == rob,],
            quantiles = q)

        MR <- quantile.table.from.samples("mr_lrt",
            PAR.SAMPLES[PAR.SAMPLES$rho==r & PAR.SAMPLES$robust == rob,],
            quantiles = q)
    
        for (n in sample_sizes) {
            nstr <- sprintf("n=%d", n)
            if (latex) {
                if (n == sample_sizes[1]) {
                    printf ("n = %d", n)
                } else {
                    printf ("%d", n)
                }
            } else {
                printf("%14s ", nstr)
            }
            print_values(RW, n)
            if (!latex) cat(" | ")
            print_values(MR, n)
    
            if (latex) {
                printf (" \\\\")
            }
            printf("\n")
        }      
    }

    qstring <- sapply(q, function(v) sprintf("p=%.2f", v))
    
    if (latex) {
        printf("\\begin{tabular}{rrrr|rrr}\n")
        printf(" & \\multicolumn{3}{c}{RW Null} & \\multicolumn{3}{c}{AR(1) Null} \\\\\n")
        printf(" RW Samples")
        qstring <- sapply(q, function(v) sprintf("%.2f", v))
        qstring[1] <- sprintf("p=%.2f", q[1])
        qstring1 <- qstring
        qstring2 <- qstring
        qstring1[length(qstring1)] <- sprintf("\\multicolumn{1}{r}{%s}", qstring1[length(qstring1)])
        vprintf("&%9s ", qstring1)
        vprintf("&%9s ", qstring2)
        printf ("\\\\ \n")
        printf (" \\hline\n")
    } else {
        printf("%14s %s\n", "", ctext("Quantiles of Likelihood Ratio Function", 65))
        printf("%14s  %30s %30s\n", "", ctext("RW Null", 30), ctext("AR(1) Null", 30))
        printf("%-14s ", "RW Samples")
        vprintf("%9s ", qstring)
        printf (" | ")
        vprintf("%9s ", qstring)
        printf ("\n")
        printf("%14s %s\n", "", paste(rep("-", 65), collapse=""))
    }
    
    print.two.panels(1)

    if (latex) {
        printf ("\\\\ \n")
        printf(" AR(1) Samples \\\\ \n")
        printf(" $\\rho = %.2f$", rho)
        vprintf("&%9s ", qstring1)
        vprintf("&%9s ", qstring2)
#        printf (" $\\rho = %.2f$ & \\multicolumn{6}{\\hline}\\\\\n", rho)
        cat("\\\\ \n")
        printf (" \\hline\n")
    } else {
        printf("\n%-14s ", "AR(1) Samples")
        vprintf("%9s ", qstring)
        printf (" | ")
        vprintf("%9s ", qstring)
        printf ("\n")
        printf("%14s %s\n", "", paste(rep("-", 65), collapse=""))
    }

    print.two.panels(rho)
    
    if (latex) {
        printf ("\\\\ \n")
        printf ("$Q[\\chi^2]$ ")
    } else {
        printf ("\n%-14s ", "Q[chi^2]")
    }
    
    for (p in q) {
        if (latex) cat (" & ")
        printf("%9.2f ", -qchisq(1 - p, 2)/2)
    }
    if (!latex) printf (" | ")
    for (p in q) {
        if (latex) cat (" & ")
        printf("%9.2f ", -qchisq(1 - p, 1)/2)
    }
    
    if (latex) {
        printf ("\\\\ \n")
        printf ("\\end{tabular}\n")
    }
    
    printf("\n\n")
    
}

par.powersamples.table <- function (nrep=1000, 
    sample_size=c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), 
    rho=c(0.1,0.3,0.5,0.7,0.9), 
    sigma_M=seq(0.5,3.0,by=.5),
    sigma_R=1,
    robust=c(FALSE, TRUE),
    use.multicore=TRUE, lambda=0,
    save.file="PAR.POWER.SAMPLES") {
	# For each value of rho, sigma_M, and sigma_R generates nrep random par series 
	# For each series, performs a fit using the maximum
    # likelihood method, and uses the likelihood ratio test to assess whether
    # or not the null hypothesis of a pure random walk should be rejected
    # and whether or not the null hypothesis of a pure AR(1) series should
    # also be rejected.
    # Returns a data.frame with the values that were generated.
	
    res <- expand.grid(1:nrep, sample_size, rho, sigma_M, sigma_R, robust)
    colnames(res) <- c("rep", "sample_size", "rho", "sigma_M", "sigma_R", "robust")
    
    do_rep <- function(k) {
        r <- res[k,]
        s <- sample.likelihood_ratio.par (r$sample_size, r$rho, r$sigma_M, r$sigma_R,
            nrep = 1, use.multicore=use.multicore, robust=r$robust, seed.start=k+1e6)
            
        s$rwnull.p.lrt <- par.rw.pvalue(s$rw_lrt, s$n, s$robust)
        s$mrnull.p.lrt <- par.mr.pvalue(s$mr_lrt, s$n, s$robust, ar1test="lr")
        s$mrnull.p.kpss <- par.mr.pvalue(-s$kpss_stat, s$n, s$robust, ar1test="kpss")
        s$par.p.lrt <- par.joint.pvalue(s$rw_lrt, s$mr_lrt, s$n, s$robust, ar1test="lr")
        s$par.p.kpss <- par.joint.pvalue(s$rw_lrt, -s$kpss_stat, s$n, s$robust, ar1test="kpss")

        s
    }

    if(use.multicore) lapply <- mclapply
    nres <- do.call("rbind", lapply(1:nrow(res), do_rep))
    if (!is.na(save.file)) write.table(nres, save.file)
	nres
}

par.plot.power <- function (which.test=c("par","rw","mr"), robust=FALSE, ar1test=c("lr","kpss"),
    rho = c(0.5,0.7,0.9), sigma_M=c(0.5,1,2)) {

    rwnull.p.lrt <- mrnull.p.lrt <- mrnull.p.kpss <-
        par.p.lrt <- par.p.kpss <- n <- sigma_R <- par.lrt <-
        rho_str <- par.kpss <- mrnull.lrt <- mrnull.kpss <-
        rwnull.lrt <- sigma <- M <- . <- x <- NULL  # Make R CMD check happy
    
    which.test <- match.arg(which.test)
    ar1test <- match.arg(ar1test)

    if (!exists("PAR.POWER.SAMPLES")) {
        if (file.exists("PAR.POWER.SAMPLES")) {
            PAR.POWER.SAMPLES <<- read.table("PAR.POWER.SAMPLES", header=TRUE)
        } else {
            PAR.POWER.SAMPLES <<- par.powersamples.table()
        }
    }
    probust <- if(robust) 1 else 0
    ADT <- data.table(PAR.POWER.SAMPLES)
    ADT <- ADT[robust == probust,]
    APS <- ADT[,list(rwnull.lrt = mean(rwnull.p.lrt < 0.05), 
                     mrnull.lrt=mean(mrnull.p.lrt < 0.05), 
                     mrnull.kpss=mean(mrnull.p.kpss < 0.05), 
                     par.lrt=mean(par.p.lrt < 0.05), 
                     par.kpss=mean(par.p.kpss < 0.05)), 
                by=list(n,rho,sigma_M,sigma_R,robust)]
                
#    partp <- part.ptab[part.ptab$sigma_M == sigma_M & part.ptab$sigma_R == sigma_R,]

    rho_list = rho
    sigma_M_list = sigma_M
    APS <- APS[rho %in% rho_list & sigma_M %in% sigma_M_list,]
               
    APS$rho_str <- sprintf("%.2f", APS$rho)
    APS$sigma_M_str <- sprintf("%.2f", APS$sigma_M)

    if (which.test == "par" && ar1test == "lr") {
        p <- ggplot (APS, aes(x = n, y=par.lrt, colour=rho_str)) + 
                ggtitle ("Power of PAR Test")
    } else if (which.test == "par" && ar1test == "kpss") {
        p <- ggplot (APS, aes(x = n, y=par.kpss, colour=rho_str)) + 
                ggtitle ("Power of PAR Test (KPSS)")
    } else if (which.test == "mr" && ar1test == "lr") {
        p <- ggplot (APS, aes(x = n, y=mrnull.lrt, colour=rho_str)) + 
                ggtitle ("Power of AR(1) vs. PAR")
    } else if (which.test == "mr" && ar1test == "kpss") {
        p <- ggplot (APS, aes(x = n, y=mrnull.kpss, colour=rho_str)) + 
                ggtitle ("Power of AR(1) vs. PAR")
    } else {
        p <- ggplot (APS, aes(x = n, y=rwnull.lrt, colour=rho_str)) + 
                ggtitle ("Power of RW vs. PAR")
    }
    
    p <- p + geom_line() +
        facet_grid (sigma_M ~ ., labeller=label_bquote(sigma[M] == .(x))) +
                ylab("Power") + xlab("Sample Size") + scale_colour_discrete(expression(rho))
    p
}

par.powersamples.pvmr.table <- function (nrep=200, 
    sample_size=c(250, 500, 750, 1000), 
    rho=seq(0.05, 0.95, by=0.05), 
    pvmr=c(0.1, 0.3, 0.5, 0.7, 0.9),
    robust=c(FALSE, TRUE),
    use.multicore=TRUE, lambda=0, 
    save.file="PAR.POWER.PVMR.SAMPLES") {
    # Generates a table of likelihood ratio tests for varying values of
    # sample_size, rho and pvmr.  
    # Returns a data.frame with the values that were generated.
    
    res <- expand.grid(sample_size, rho, pvmr, robust)
    colnames(res) <- c("sample_size", "rho", "pvmr","robust")
    
    do_rep <- function(k) {
        r <- res[k,]        
        sigma_M <- sqrt((1 + r$rho) * r$pvmr / (2 - 2 * r$pvmr))
        sigma_R <- 1
        df <- par.powersamples.table (nrep, r$sample_size, r$rho, 
            sigma_M, sigma_R, robust=r$robust)        
        df$pvmr <- r$pvmr
        df
    }

    ndf <- do.call("rbind", lapply(1:nrow(res), do_rep))
    if (!is.na(save.file)) write.table(ndf, save.file)
	ndf   
}

par.plot.power.pvmr <- function (robust=FALSE, ar1test=c("lr","kpss"),
    pvmr=c(0.3, 0.5, 0.7)) {

    rwnull.p.lrt <- mrnull.p.lrt <- mrnull.p.kpss <-
        par.p.lrt <- par.p.kpss <- n <- rho <- sigma_M <- sigma_R <-
        par.lrt <- R2 <- par.kpss <- . <- x <- NULL # Make R CMD check happy
    
#    partp <- part.ptab[part.ptab$sigma_M == sigma_M & part.ptab$sigma_R == sigma_R,]
#    partp$rho <- sprintf("%.2f", partp$rho)
#    partp$sigma_M <- sprintf("%.2f", partp$sigma_M)

    ar1test <- match.arg(ar1test)
    if (!exists("PAR.POWER.PVMR.SAMPLES")) {
        if (file.exists("PAR.POWER.PVMR.SAMPLES")) {
            PAR.POWER.SAMPLES <<- read.table("PAR.POWER.PVMR.SAMPLES", header=TRUE)
        } else {
            PAR.POWER.SAMPLES <<- par.powersamples.pvmr.table()
        }
    }
    probust <- if(robust) 1 else 0
    ADT <- data.table(PAR.POWER.PVMR.SAMPLES)
    ADT <- ADT[robust == probust,]
    pvmr.ptab <- ADT[,list(rwnull.lrt = mean(rwnull.p.lrt < 0.05), 
                     mrnull.lrt=mean(mrnull.p.lrt < 0.05), 
                     mrnull.kpss=mean(mrnull.p.kpss < 0.05), 
                     par.lrt=mean(par.p.lrt < 0.05), 
                     par.kpss=mean(par.p.kpss < 0.05),
                     pvmr=mean(pvmr)), 
                by=list(n,rho,sigma_M,sigma_R,robust)]
    
    pvmr.ptab$R2 <- sprintf("%.2f", pvmr.ptab$pvmr)
    pvmr.ptab$sample_size <- sprintf("%4d", pvmr.ptab$n)

    pvmr_list <- pvmr
    pvmr.ptab <- pvmr.ptab[pvmr %in% pvmr_list,]
               
    pvmr.ptab$rho_str <- sprintf("%.2f", pvmr.ptab$rho)
    pvmr.ptab$sigma_M_str <- sprintf("%.2f", pvmr.ptab$sigma_M)

    if (ar1test == "lr") {
        p <- ggplot (pvmr.ptab, aes(x = rho, y=par.lrt, colour=R2))
    } else {
        p <- ggplot (pvmr.ptab, aes(x = rho, y=par.kpss, colour=R2))
    }
            
    
    p <- p + geom_line() +
        facet_grid (sample_size ~ ., labeller=label_bquote("Sample Size" == .(x))) +
                ylab("Power") + xlab(expression(rho)) +
                scale_colour_discrete(name = bquote(R[MR]^2)) +
                ggtitle ("Power of PAR Test")
    p
}

par.plot.lrt.vs.kpss <- function (robust=FALSE, 
    rho = c(0.5,0.7,0.9), sigma_M=c(0.5,1,2)) {

    rwnull.p.lrt <- NULL # Make R CMD check happy
    mrnull.p.lrt <- NULL # Make R CMD check happy
    par.p.lrt <- NULL # Make R CMD check happy
    par.p.kpss <- NULL # Make R CMD check happy
    n <- NULL # Make R CMD check happy
    sigma_R <- NULL # Make R CMD check happy
    label <- NULL # Make R CMD check happy
    mrnull.p.kpss <- NULL
    
    if (!exists("PAR.POWER.SAMPLES")) {
        PAR.POWER.SAMPLES <<- read.table("PAR.POWER.SAMPLES", header=TRUE)
    }
    probust <- if(robust) 1 else 0
    ADT <- data.table(PAR.POWER.SAMPLES)
    ADT <- ADT[robust == probust,]
    APS <- ADT[,list(rwnull.lrt = mean(rwnull.p.lrt < 0.05), 
                     mrnull.lrt=mean(mrnull.p.lrt < 0.05), 
                     mrnull.kpss=mean(mrnull.p.kpss < 0.05), 
                     par.lrt=mean(par.p.lrt < 0.05), 
                     par.kpss=mean(par.p.kpss < 0.05)), 
                by=list(n,rho,sigma_M,sigma_R,robust)]
            
#    partp <- part.ptab[part.ptab$sigma_M == sigma_M & part.ptab$sigma_R == sigma_R,]

    rho_list = rho
    sigma_M_list = sigma_M
    APS <- APS[rho %in% rho_list & sigma_M %in% sigma_M_list,]
    APS$rho_str <- sprintf("rho == %.2f", APS$rho)
    APS$sigma_M_str <- sprintf("sigma[M] == %.2f", APS$sigma_M)

    df1 <- APS
    df1$value <- df1$par.lrt
    df1$label <- "LR"
    
    df2 <- APS
    df2$value <- df2$par.kpss
    df2$label <- "KPSS"
    
    df <- rbind(df1, df2)
    p <- ggplot(df, aes(x = n, y=value, colour=label)) +
        geom_line() + ylab("Power") + xlab("Sample Size") +
        facet_grid (rho_str ~ sigma_M_str, labeller=label_parsed) +
        ggtitle("Comparison of Likelihood Ratio (LR) to KPSS Test")
    
    p
}
