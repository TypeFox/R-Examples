#################################
### Main estimation functions:
#################################

dod <-
    function(same, diff, d.prime0 = 0, conf.level = 0.95,
             statistic = c("likelihood", "Pearson", "Wilcoxon", "Wald"),
             alternative = c("difference", "similarity", "two.sided",
             "less", "greater"), control=dodControl(), ...)
### Our main fitting function.
### ... parsed on to wilcox.test when statistic = "Wilcoxon".
{
    ## Match arg:
    alt <- match.arg(alternative)
    if(alt == "difference") alt <- "greater"
    if(alt == "similarity") alt <- "less"
    stat <- match.arg(statistic)
    ctrl <- control
    if(!inherits(control, "dodControl"))
        stop("Specify 'control' with dodControl()")
    ## Test arg:
    if(ctrl$test.args)
        testdodargs(same=NULL, diff=NULL, conf.level=conf.level,
                    d.prime0=d.prime0,
                    return.data=FALSE, int.tol=ctrl$integer.tol)
    if(isTRUE(all.equal(d.prime0, 0)) &&
       !alternative %in% c("difference", "greater"))
        stop("'alternative' has to be 'difference' or 'greater' if 'd.prime0' is 0")
    if(!isTRUE(all.equal(d.prime0, 0)) && stat == "Wilcoxon") {
        stop("Wilcoxon statistic only available with d.prime0 = 0")
    }
    ## Fit DOD model:
    fit <- dod_fit(same, diff, control=ctrl, ...)
    ## Format data:
    .data <- fit$data
    nlev <- ncol(.data)
    npar <- length(fit$coefficients)
    same <- as.vector(.data[1, ])
    diff <- as.vector(.data[2, ])
    logLik <- fit$logLik
    d.prime <- fit$d.prime
    ## Compute standard errors:
    if(d.prime < 0.01 && ctrl$get.vcov && ctrl$do.warn) {
        warning("d.prime < 0.01: standard errors are unavailable")
    }
    std.err <- rep(NA_real_, nlev)
    if(!is.null(fit$vcov) && all(is.finite(fit$vcov)))
        std.err <- sqrt(diag(fit$vcov))
    names(std.err) <- names(fit$coefficients)
    ## Compute confidence intervals:
    if(stat == "Wald") {
        ci <- confint(fit, level=conf.level, method="Wald")
        conf.method <- "Wald"
    } else {
        ci <- confint(fit, level=conf.level, method="likelihood")
        conf.method <- "profile likelihood"
    }
### FIXME: add convergence information to model, so it can be
### accessed and assessed.
    ## Collect estimates in tables:
    coef.tab <- cbind(d.prime, std.err[npar], ci[1L], ci[2L])
    rownames(coef.tab) <- "d.prime"
    colnames(coef.tab) <- c("Estimates", "Std. Error", "Lower", "Upper")
    boundary.coeftab <- rbind(fit$tau, std.err[-npar])
    rownames(boundary.coeftab) <- c("Estimate", "Std. Error")
    colnames(boundary.coeftab) <-
        as.character(seq_len(ncol(boundary.coeftab)))
    ## Results list:
    fit$coefficients <- coef.tab
    res <- c(fit,
             list(conf.level=conf.level, boundary.coeftab=boundary.coeftab,
                  statistic=stat, d.prime0=d.prime0, alternative=alt,
                  conf.method=conf.method))
    res$call <- match.call(expand.dots=TRUE)
    ## Get p-value and value of test statistic:
    res$stat.value <- res$p.value <- NA_real_
    if(stat == "Wilcoxon") {
        ## Assumes that all response categories have at least one
        ## positive record:
        wtest <- wilcox.test(rep.int(1:nlev, diff), rep.int(1:nlev, same),
                             alternative=alt, ...)
        res$stat.value <- unname(wtest$statistic)
        res$p.value <- unname(wtest$p.value)
    } else if(stat == "Wald") {
        if(is.finite(std.err[npar])) {
            res$stat.value <- unname((d.prime - d.prime0) / std.err[npar])
            res$p.value <- unname(normalPvalue(res$stat.value, alternative=alt))
        }
    } else if(stat %in% c("likelihood", "Pearson")) {
        ## Fit model under the null hypothesis where d.prime = d.prime0:
        ctrl$get.vcov <- ctrl$get.grad <- ctrl$test.args <- FALSE
        fit0 <- dod_fit(same, diff, d.prime=d.prime0, control=ctrl, ...)
        logLik0 <- fit0$logLik
        d.prime0 <- fit0$d.prime
        ## Assess convergence of model under the alternative:
        if(logLik < logLik0 && abs(logLik - logLik0) > 1e-6 && ctrl$do.warn)
            warning("Estimation of DOD model failed: likelihood did not increase")
        if(stat == "likelihood") {
            LR <- 2 * (logLik - logLik0) ## LR shold be positive
            LR[LR < 0 && abs(LR) < 1e-4] <- 0
            if(LR >= 0) {
                res$stat.value <- unname(sign(d.prime - d.prime0) * sqrt(LR))
                res$p.value <- unname(normalPvalue(res$stat.value, alternative=alt))
            }
        } else if(stat == "Pearson") {
            ## Pearson equivalent of LR test:
            freq <- par2prob_dod(tau=fit$tau, d.prime=d.prime) * rowSums(.data)
            ## Probabilities and frequencies under the null model (d.prime0):
            freq0 <- par2prob_dod(tau=fit0$tau, d.prime=fit0$d.prime) *
                rowSums(.data)
### FIXME: Check that this is the right test statistic (or should 'freq' be in the
### denominator?)
            X2 <- sum((freq - freq0)^2 / freq0)
            res$stat.value <- unname(sign(d.prime - d.prime0) * sqrt(X2))
            res$p.value <- unname(normalPvalue(res$stat.value, alternative=alt))
        }
    } else stop("'statistic' not recognized")
    ## Return:
    class(res) <- "dod"
    res
}

dod_fit <-
    function(same, diff, tau=NULL, d.prime=NULL, control=dodControl(),
             ...)
{
    ## Return list:
    res <- list(call = match.call(expand.dots=TRUE))
    res$control <- control
    if(!inherits(control, "dodControl"))
        stop("Specify 'control' with dodControl()")
    ## Read data:
    x <- readdoddata(same, diff)
    ## Test arg:
    if(control$test.args)
        x <- testdodargs(x$same, x$diff, tau=tau, d.prime0=d.prime,
                         int.tol=res$control$integer.tol,
                         return.data=TRUE)
    same <- x$same
    diff <- x$diff
    nlev <- length(same)
    .data <- rbind(same, diff)
    colnames(.data) <- as.character(seq_len(ncol(.data)))
    rownames(.data) <- c("same-pairs", "diff-pairs")
    res$data <- .data
    res$tau <- as.vector(tau)
    res$d.prime <- as.vector(d.prime)
### 3 CASES to get parameters:
    if(is.null(tau) && !is.null(d.prime) &&
       isTRUE(all.equal(d.prime, 0))) {
### CASE 1:
        ## d-prime = 0 and tau=NULL:
        ## No optimization necessary.
        res$tau <- dod_null_tau_internal(same=same, diff=diff)
        res$d.prime <- d.prime
    } else if(is.null(tau) || is.null(d.prime)) {
### CASE 2:
        ## We need to optimize either tau, d.prime or c(tau, d.prime):
        if(is.null(d.prime) && is.null(tau)) { ## get c(tau, d.prime)
            start <- c(cumsum(init_tau(nlev)), 1)
            par.names <-
                c(paste("tau", seq_len(nlev-1L), sep=""), "d.prime")
            npar <- length(start)
            objfun <- function(par) {
                dod_nll_internal(tau=par[1:(npar - 1L)], d.prime=par[npar],
                                 same=same, diff=diff)
            }
            case <- "21"
        } else if(is.null(tau) && !is.null(d.prime)) { ## get tau
            start <- cumsum(init_tau(nlev))
            par.names <-  paste("tau", seq_len(nlev-1L), sep="")
            objfun <- function(par)
                dod_nll_internal(tau=par, d.prime=d.prime,
                                 same=same, diff=diff)
            case <- "22"
        } else if(!is.null(tau) && is.null(d.prime)) { ## get d.prime
            start <- 1
            par.names <- "d.prime"
            objfun <- function(par)
                dod_nll_internal(tau=tau, d.prime=par,
                                 same=same, diff=diff)
            case <- "23"
        }
        res$start <- setNames(start, par.names)
        npar <- length(start) ## Perhaps not needed
        ## Optimization:
        res$optRes <- nlminb(res$start, objfun, lower=0, control=control$optCtrl)
        par <- as.vector(res$optRes$par)
        ## Convergence tests:
        start.value <- objfun(res$start)
        if(all(par > 1e-2) && control$get.grad){
            g <- grad(objfun,  x=par)
            res$gradient <- setNames(g, par.names)
            if(!all(is.finite(g)) && control$do.warn) {
                warning("Cannot assess convergence: non-finite gradient")
            } else {
                if(max(abs(g)) > control$grad.tol && control$do.warn)
                    warning("Estimation failed with max(gradient) = ",
                            format(max(abs(g)), digits=2), " (grad.tol = ",
                            format(control$grad.tol, digits=2), ")")
                if(control$get.vcov) {
                    ## Hessian and vcov:
                    h <- res$vcov <- array(NA_real_, dim=c(npar, npar),
                                      dimnames=list(par.names, par.names))
                    h <- try(hessian(objfun, x=par), silent=TRUE)
                    if((inherits(h, "try-error") || !all(is.finite(h)))
                       && control$do.warn) {
                        warning("unable to compute Hessian")
                    } else {
                        ch <- try(chol(h), silent=TRUE)
                        if(inherits(ch, "try-error") && control$do.warn) {
                            warning("Model is ill-defined and may not have converged")
                        } else {
                            res$vcov[] <- chol2inv(ch)
                        }
                        dimnames(h) <- list(par.names, par.names)
                        res$Hessian <- h
                    }
                }
            }
        }
        ## Match parameters:
        if(case == "21") {
            res$tau <- par[1:(npar-1L)]
            res$d.prime <- par[npar]
        } else if(case == "22") {
            res$tau <- par
            res$d.prime <- as.vector(d.prime)
        } else { ## case == "23"
            res$tau <- as.vector(tau)
            res$d.prime <- par
        }
    }
### CASE 3:
    ## Both tau and d.prime set, so we just evaluate the likelihood
    ## No optimization necessary
    res$logLik <- -dod_nll_internal(tau=res$tau, d.prime=res$d.prime,
                                    same=same, diff=diff)
    ## Set parameter names:
    res$coefficients <-
        setNames(c(res$tau, res$d.prime),
                 c(paste("tau", seq_along(res$tau), sep=""), "d.prime"))
    ## Return fit:
    class(res) <- "dod_fit"
    res
}

#################################
### dodControl:
#################################

dodControl <- function(grad.tol = 1e-4,
                       ## eigen.value.tol=1e-4,
                       integer.tol = 1e-8,
                       get.vcov = TRUE,
                       get.grad = TRUE,
                       test.args = TRUE,
                       do.warn=TRUE,
                       optCtrl=list())
{
    stopifnot(length(grad.tol) == 1L,
              grad.tol > 0,
              is.finite(grad.tol))
    stopifnot(length(integer.tol) == 1L,
              integer.tol > 0,
              is.finite(integer.tol))
    stopifnot(length(get.vcov) == 1L,
              is.logical(get.vcov))
    stopifnot(length(get.grad) == 1L,
              is.logical(get.grad))
    stopifnot(length(test.args) == 1L,
              is.logical(test.args))
    stopifnot(length(do.warn) == 1L,
              is.logical(do.warn))
    ## Return value:
    res <- list(grad.tol=grad.tol,
                integer.tol=integer.tol,
                get.grad=get.grad,
                get.vcov=get.vcov,
                test.args=test.args,
                do.warn=do.warn,
                optCtrl=optCtrl)
    class(res) <- "dodControl"
    res
}

#################################
### Optimal tau estimation:
#################################

optimal_tau <-
    function(d.prime, d.prime0 = 0, ncat=3,
             method=c("equi.prob", "LR.max", "se.min"),
             tau.start=NULL, equi.tol = 1e-4, grad.tol = 1e-2,
             do.warn=TRUE)
### FIXME: Could add "Pearson.max" as an option here.
### Estimate optimal boundary parameters, tau in a DOD model with one
### of the following optimality criteria:
### 1) Equal category probabilities averaged over same-pairs and diff-pairs.
### 2) Minimum standard error of d.prime/Maximum Wald test.
### 3) Maximum LR statistics for the test: H0: d.prime = 0, HA: d.prime > 0.
{
    ## Define objective functions:
    nlr_stat <- function(tpar) {
        ## Negative LR statistic.
        tau <- cumsum(tpar)
        ## Limiting distribution of data given c(tau, d.prime):
        data <- par2prob_dod(tau=tau, d.prime=d.prime) * 100
        ## logLiks at d.prime=d.prime and d.prime=0:
        logLik <- -dod_nll_internal(tau=tau, d.prime=d.prime, same=data[1, ], diff=data[2, ])
        if(isTRUE(all.equal(d.prime0, 0)))
            logLik0 <- -dod_null_internal(same=data[1, ], diff=data[2, ])
        else
            logLik0 <- dod_fit(same=data[1, ], diff=data[2, ],
                               d.prime=d.prime0, control=ctrl)$logLik
        -(logLik - logLik0) ## negative LR/2 statistic
    }
    min_se <- function(tpar) {
        ## function to minimize se(d.prime):
        if(any(!is.finite(tpar)))
            return(Inf)
        tau <- cumsum(tpar)
        data <- par2prob_dod(tau, d.prime) * 100
        h <- try(hessian(dod_nll_all_internal, c(tau, d.prime), same=data[1, ],
                         diff=data[2, ]), silent=TRUE)
        if(inherits(h, "try-error") || any(!is.finite(h)))
            return(Inf)
        vcov <- try(solve(h), silent=TRUE)
        if(inherits(vcov, "try-error") || any(!is.finite(vcov)))
            return(Inf)
         sqrt(rev(diag(solve(h)))[1])
    }
    equi_prob <- function(tpar, const=1) {
        tau <- cumsum(tpar)
        prob <- colSums(par2prob_dod(tau, d.prime)) / 2
        const + sum(((prob  - 1/ncat) * 1e3)^2)
    }
    tau2tpar <- function(tau) c(tau[1], diff(tau))
    ## tpar2tau <- function(tpar) cumsum(tpar)
    ## Match and test arguments:
    method <- match.arg(method)
    stopifnot(is.numeric(d.prime),
              length(d.prime) == 1L,
              d.prime >= 0) ## OBS: strictly larger than 0!
    stopifnot(is.numeric(ncat),
              length(ncat) == 1L,
              ncat >= 2)
    stopifnot(length(do.warn) == 1L,
              is.logical(do.warn))
    if(abs(ncat - round(ncat)) > 1e-8 && do.warn) {
        warning("Non-integer ncat rounded to ", round(ncat))
    }
    ncat <- as.integer(round(ncat))
    ## Match objective function:
    objfun <- switch(method, "equi.prob" = equi_prob,
                     "se.min" = min_se,
                     "LR.max" = nlr_stat)
    if(method == "LR.max") {
        ctrl <- dodControl(test.args=FALSE, do.warn=FALSE,
                           get.grad=FALSE, get.vcov=FALSE,
                           integer.tol=10,
                           grad.tol=grad.tol)
    }
    ## Get starting values:
    if(!is.null(tau.start)) { ## 1, 2, 3 -> c(ts[1], diff(ts))
        stopifnot(is.numeric(tau.start),
                  all(is.finite(tau.start)),
                  all(tau.start > 0))
        if(any(diff(tau.start) <= 0))
            stop("'tau.start' is not increasing")
        if(length(tau.start) != ncat-1L)
            stop("length(tau.start) is ", length(tau.start),
                 ": was expecting length(tau.start) = ncat - 1 = ",
                 ncat - 1L)
        start <- tau2tpar(tau.start)
    } else
        start <- init_tau(ncat)
    ## Estimate optimal tau:
    fit <- nlminb(start=start, objective=objfun, lower=1e-4)
    ## if(fit$convergence != 0)
    ##     warning("Estimation of tau failed to converge with nlminb code ",
    ##             fit$convergence)
    tau <- cumsum(fit$par)
    prob <- par2prob_dod(tau, d.prime)
    if(method == "equi.prob") {
        diffs <- abs((colSums(prob) / 2) - 1/ncat)
        if(max(diffs) > equi.tol && do.warn)
            warning("Estimation of tau failed with max(diffs) = ",
                    format(max(diffs), digits=2), " (equi.tol = ",
                    format(equi.tol, digits=2), ")")
    }
    ## Obtain and check gradient:
    g <- grad(func=objfun, x=fit$par)
    if(max(abs(g)) > grad.tol && do.warn)
        warning("Estimation of tau failed with max(gradient) = ",
                format(max(abs(g)), digits=2), " (grad.tol = ",
                format(grad.tol, digits=2), ")")
    ## Results:
    res <- list(tau=tau, prob=prob, tau.start=cumsum(start), gradient=g,
                fit=fit, method=method, call=match.call(expand.dots=TRUE),
                grad.tol=grad.tol, equi.tol=equi.tol)
    res
}

#################################
### Misc estimation functions:
#################################

dod_null_tau_internal <- function(same, diff) {
### Compute tau-vector for the DOD null model (d.prime=0)
    alldata <- rev(same + diff)
    nlev <- length(alldata)
    cumdata <- cumsum(alldata)
    tau <- sapply(1:(nlev - 1), function(lev) {
        -sqrt(2) * qnorm(cumdata[lev] / (2* sum(alldata)))
    })
    tau <- rev(tau)
    tau
}

dod_null_tau <- function(same, diff) {
### Compute tau-vector for the DOD null model (d.prime=0)
    x <- readdoddata(same, diff)
    x <- testdodargs(x$same, x$diff, int.tol=1) ## don't warn with non-integers
    dod_null_tau_internal(x$same, x$diff)
}

par2prob_dod_internal <- function(tau, d.prime) {
    ## Cumulative probilities:
    gamma.same <- 2*pnorm(tau / sqrt(2)) - 1
    gamma.diff <- pnorm((tau - d.prime)/sqrt(2)) -
        pnorm((-tau - d.prime)/sqrt(2))
    ## Non-cumulative probabilities:
    p.same <- c(gamma.same, 1) - c(0, gamma.same)
    p.diff <- c(gamma.diff, 1) - c(0, gamma.diff)
    ## Make sure sum(p.same) == sum(p.diff) == 1:
    p.same <- p.same / sum(p.same)
    p.diff <- p.diff / sum(p.diff)
    ## Return prob matrix:
    rbind(p.same, p.diff)
}

par2prob_dod <- function(tau, d.prime) {
### Computes the dod probability matrix from dod parameters tau and
### d.prime; rbind(same-pair probs, diff-pair probs);
### Result satisfies: 1) rowSums(<return_value>) == c(1, 1),
### 2) ncol(<return_value>) == length(tau) + 1L.
    stopifnot(length(d.prime) == 1L,
              is.numeric(d.prime),
              d.prime >= 0)
    stopifnot(is.numeric(tau),
              length(tau) >= 1L,
              all(tau > 0),
              all(diff(tau) > 0))
    par2prob_dod_internal(tau, d.prime)
}

init_tau <- function(ncat=4) c(1, rep(3/(ncat-1), ncat-2))

#################################
### Likelihood functions:
#################################

dod_null_internal <- function(same, diff) {
### Compute neg. logLik for DOD null model (with d.prime=0)
    tau <- dod_null_tau(same, diff)
    dod_nll_internal(tau=tau, d.prime=0, same=same, diff=diff)
}

dod_null <- function(same, diff, integer.tol=1e-8) {
### Like dod_null_internal, but also tests admissibility of
### arguments.
### Computes the neg. logLik for the DOD null model (with d.prime=0).
### same and diff can be factors or numeric vectors.
    ## Read same & diff:
    x <- readdoddata(same, diff)
    ## Test same and diff args:
    x <- testdodargs(x$same, x$diff, int.tol=integer.tol)
    ## Compute logLik:
    dod_null_internal(x$same, x$diff)
}

dod_nll_internal <- function(tau, d.prime, same, diff) {
### Computes the nll of the dod Thurstonian model.
### Returns Inf if parameters are inadmissible.
    prob <- par2prob_dod_internal(tau, d.prime)
    ## compute nll:
    if(all(is.finite(prob)) && all(prob > 0))
    ## if(all(prob > 0))
        -sum(same %*% log(prob[1, ])) -
            sum(diff %*% log(prob[2, ]))
    else
        Inf
}

dod_nll <- function(tau, d.prime, same, diff, integer.tol=1e-8) {
### Computes the neg. logLik of the dod Thurstonian model.
### Returns an error if the parameters are inadmissible.

    ## Get case probabilities:
    prob <- par2prob_dod(tau, d.prime)
### NOTE: tau and d.prime args are tested in par2prob_dod().
    if(any(prob <= 0)) stop("Inadmissible parameters")
    ## Read same & diff:
    x <- readdoddata(same, diff)
    ## Test same and diff args:
    x <- testdodargs(x$same, x$diff, tau=tau, int.tol=integer.tol)
    ## Compute negative log-likelihood:
    nll <- -sum(x$same %*% log(prob[1, ])) -
        sum(x$diff %*% log(prob[2, ]))
    nll
}

dod_nll_all_internal<- function(par, same, diff) {
### nll expressed as a function of all=c(tau, d.prime) parameters.
    npar <- length(par)
    dod_nll_internal(tau=par[1:(npar-1)], d.prime=par[npar], same, diff)
}

#################################
### Utils:
#################################

readdoddata <- function(same, diff) {
### Read same & diff arguments to dod() etc.
### Returns a list of same and diff numeric vectors of observations.
    if(is.factor(same) && is.factor(diff)) {
        stopifnot(nlevels(same) == nlevels(diff))
        stopifnot(all(levels(same) == levels(diff)))
        same <- as.vector(table(same))
        diff <- as.vector(table(diff))
    }
    else {
        same <- as.vector(same)
        diff <- as.vector(diff)
    }
    list(same=same, diff=diff)
}

testdodargs <-
    function(same=NULL, diff=NULL, tau=NULL, conf.level=NULL,
             d.prime0=NULL, int.tol=1e-8, return.data=TRUE)
{
### Tests admissibility of same and diff arguments.
    stopifnot(length(return.data) == 1L,
              is.logical(return.data))
    ## Test integer.tol argument:
    stopifnot(is.numeric(int.tol),
              length(int.tol) == 1L,
              int.tol >= 0)
    ## Test same & diff:
    has.data <- !is.null(same) && !is.null(diff)
    if(has.data) {
        stopifnot(is.numeric(same), is.numeric(diff),
                  all(same >= 0),
                  all(diff >= 0),
                  length(same) == length(diff))
        if(!is.null(tau))
            stopifnot(length(same) == length(tau) + 1L)
        .data <- rbind(same, diff)
        .data <- .data[, colSums(.data) > 0, drop=FALSE]
### FIXME: warn if any(colSums(.data) <= 0)?
        if(ncol(.data) <= 1L)
            stop("Need counts in more than one response category", call.=FALSE)
        if(any(abs(.data - round(.data)) > int.tol))
            warning("non-integer counts in 'same' or 'diff'", call.=FALSE)
    }
    ## conf.level:
    if(!is.null(conf.level)) {
        stopifnot(is.numeric(conf.level),
                  length(conf.level) == 1,
                  0 < conf.level,
                  conf.level < 1)
    }
    ## d.prime0:
    if(!is.null(d.prime0)) {
        stopifnot(is.numeric(d.prime0),
                  length(d.prime0) == 1,
                  d.prime0 >= 0)
    }
    ## Return modified data:
    if(has.data && return.data)
        list(same=as.vector(drop(.data[1, ])),
             diff=as.vector(drop(.data[2, ])))
    else
        invisible()
}

#################################
### Print methods:
#################################

print.dod <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat(paste("\nResults for the Thurstonian model for the Degree-of-Difference method",
              "\n\n"))
    cat(paste("Confidence level for 2-sided ", x$conf.method,
              " interval: ", round(100 * x$conf.level, 3), "%\n\n",
              sep=""))
    print(x$coefficients, digits=digits)
    cat("\nBoundary coefficients:\n")
    print(x$boundary.coeftab, digits=digits)
    cat("\nData:\n")
    print(as.data.frame(x$data))

    cat("\nResults of discrimination test:\n")
    if(x$statistic == "Wald")
        cat(paste("Wald statistic = ", format(x$stat.value, digits),
                  ", p-value = ", format.pval(x$p.value, digits=4), "\n",
                  sep=""))
    if(x$statistic == "likelihood")
        cat(paste("Likelihood Root statistic = ",
                  format(x$stat.value, digits),
                  ", p-value = ", format.pval(x$p.value, digits=4), "\n",
                  sep=""))
    if(x$statistic == "Pearson")
        cat(paste("Signed Pearson root statistic = ",
                  format(x$stat.value, digits),
                  ", p-value = ", format.pval(x$p.value, digits=4), "\n",
                  sep=""))
    if(x$statistic == "Wilcoxon")
        cat(paste("Wilcoxon statistic = ", x$stat.value,
                  "p-value =", format.pval(x$p.value, digits=4), "\n"))
    cat("Alternative hypothesis: ")
    alt <- x$alternative
    cat(paste("d-prime is",
              ifelse(alt == "two.sided", "different from",
                     paste(alt, "than")), format(x$d.prime0, digits),
              "\n"))
    invisible(x)
}

#################################
### Simulation functions:
#################################

get_tau <-
    function(d.prime, tau=NULL, ncat=NULL,
             method.tau = c("LR.max", "equi.prob", "se.min", "user.defined"),
             d.prime0=0, ...)
### Internal function.
### Computes tau if not already present and tests that arguments
### conform.
{
    method <- match.arg(method.tau)
    ## Test conformance of tau and method.tau:
    if(method != "user.defined") {
        if(!is.null(tau))
            warning("'tau' is ignored when method.tau != 'user.defined' (method.tau was '",
                    method, "')")
        stopifnot(length(ncat) == 1L,
                  is.numeric(ncat),
                  ncat >= 2L,
                  is.finite(ncat))
        if(abs(ncat - round(ncat)) > 1e-6)
            warning("non-integer 'ncat': ", format(ncat, digits=8),
                    " rounded to ", round(ncat))
        ncat <- round(ncat)
        tau <- optimal_tau(d.prime=d.prime, ncat=ncat, method=method,
                           do.warn=FALSE, d.prime0=d.prime0, ...)$tau ## grad.tol, equi.tol
    } else { ## method == "user.defined"
        if(is.null(tau))
            stop("specify 'tau' when method.tau == 'user.defined'")
        stopifnot(length(tau) >= 1L,
                  is.numeric(tau),
                  all(diff(tau) > 0),
                  all(is.finite(tau)))
### FIXME: Do we need more tests of tau here?
    }
    ## Return:
    as.vector(tau)
}

dodSim <-
    function(d.prime, ncat=4, sample.size = c(100, 100),
             method.tau = c("equi.prob", "LR.max", "se.min", "user.defined"),## max.X2?
             tau = NULL, d.prime0 = 0, ...)
{
    ## Test d.prime:
    stopifnot(length(d.prime) == 1L,
              is.numeric(d.prime),
              d.prime >= 0,
              is.finite(d.prime))
    stopifnot(length(d.prime0) == 1L,
              is.numeric(d.prime0),
              d.prime0 >= 0,
              is.finite(d.prime0))
    ## Match and test tau-related args:
    method <- match.arg(method.tau)
    tau <- get_tau(d.prime=d.prime, tau=tau, ncat=ncat,
                   method.tau=method, d.prime0=d.prime0)

    ## Test sample.size arg:
    stopifnot(is.numeric(sample.size),
              all(sample.size >= 0),
              length(sample.size) == 1L || length(sample.size) == 2L)
### NOTE: sample.size == 0 is allowed and runs with expected results.
    size <- sample.size
    if(max(abs(size - round(size))) > 1e-3)
        warning("Non-integer 'sample.size' rounded to ", round(size))
    size <- round(size)
    if(length(size) == 1L) size <- c(size, size)

    ## Get prob and simulate data:
    prob <- par2prob_dod(tau, d.prime)
### NOTE: tau and d.prime args are tested in par2prob_dod().
    mat <- t(cbind(rmultinom(n=1, size=size[1], prob=prob[1, ]),
                   rmultinom(n=1, size=size[2], prob=prob[2, ])))
    dimnames(mat) <- list(c("same-pairs", "diff-pairs"),
                          seq_len(ncol(mat)))

    ## Return:
    mat
}


#################################
### profile likelihood CI
#################################

lroot <- function(d.prime, same, diff) {
    ## Profile likelihood root statistic for d.prime
    ## Currently not used.
    ctrl <- dodControl(get.grad=FALSE, get.vcov=FALSE, do.warn=FALSE,
                       test.args=FALSE)
    fit <- dod_fit(same, diff, control=ctrl)
    logLik <- fit$logLik
    mle <- fit$d.prime
    logLik0 <- dod_fit(same, diff, d.prime=d.prime, control=ctrl)$logLik
    lroot <- sign(mle - d.prime) * sqrt(2 * abs(logLik - logLik0))
    lroot
}

get_ci <-
    function(same, diff, alpha = 0.05, dod.fit=NULL,
             tol.d.prime = 1e-3)
{
### Computes the profile likelihood CI for d-prime.
    lroot_fun <- function(d.prime, lim) {
        ## Profile likelihood root statistic minus lim.
        ll <- dod_fit(same=same, diff=diff, d.prime=d.prime,
                      control=ctrl)$logLik
        lroot <- sign(mle - d.prime) * sqrt(2 * abs(dod.fit$logLik - ll))
        lroot - lim
    }

    lim <- qnorm(p=1-alpha/2, lower.tail=TRUE)
    ctrl <- dodControl(get.vcov=FALSE, get.grad=FALSE, do.warn=FALSE,
                       test.args=FALSE)
    ## Get fit if missing:
    if(is.null(dod.fit)) {
        dod.fit <- dod_fit(same, diff, control=ctrl)
    } else {
        stopifnot(class(dod.fit) %in% c("dod", "dod_fit"))
    }
    mle <- dod.fit$d.prime
    ## Get lower and upper limits for the lroot search:
    lwr <- lroot_lwr(same, diff, mle=mle, lim=lim,
                     logLik=dod.fit$logLik, tol.d.prime=tol.d.prime)
    if(length(lwr) == 2)
        lwr <- uniroot(f=lroot_fun, interval=lwr, tol=1e-4, lim=lim)$root
    upr <- lroot_upr(same, diff, mle=mle, lim=lim,
                     logLik=dod.fit$logLik, tol.d.prime=tol.d.prime)
    if(length(upr) == 2)
        upr <- uniroot(f=lroot_fun, interval=upr, tol=1e-4, lim=-lim)$root
    ## Return CI:
    ci <- as.vector(c(lwr, upr))
    ci
}

lroot_upr <- function(same, diff, mle, lim, logLik, tol.d.prime=1e-3, ...)
{
    if(is.na(mle)) return(NA_real_)
    if(mle == Inf) return(Inf)
    ctrl <- dodControl(get.vcov=FALSE, get.grad=FALSE, do.warn=FALSE,
                       test.args=FALSE)
    Upr <- 10
    if(mle > Upr) {
        warning("Cannot determine upper confidence bound")
        return(NA_real_)
    } else {
        ll.Upr <- dod_fit(same=same, diff=diff, d.prime=Upr,
                          control=ctrl)$logLik
        LR.Upr <- 2 * (logLik - ll.Upr)
        if(LR.Upr < lim^2) {
            warning("Cannot determine upper confidence bound")
            return(NA_real_)
        }
    }
    c(mle, Upr)
}

lroot_lwr <- function(same, diff, mle, lim, logLik, tol.d.prime=1e-3, ...)
{
    ctrl <- dodControl(test.args=FALSE, get.vcov=FALSE,
                       get.grad=FALSE, do.warn=FALSE)
    if(is.na(mle)) return(NA_real_)
    if(mle < tol.d.prime) return(0)

    LR0 <- 2 * (logLik + dod_null_internal(same, diff))
    if(LR0 < lim^2) return(0)
    ## positive lwr exists:
    Upr <- mle
    if(mle == Inf) {
        ## determine upper bound for search
        LR10 <- 2 * (logLik +
                     dod_fit(same=same, diff=diff, d.prime=10,
                             control=ctrl)$logLik)
        if(LR10 > lim^2) {
            warning("Cannot determine lower confindence bound")
            return(NA_real_)
        } else {
            Upr <- 10
        }
    }
    c(0, Upr)
}

confint.dod_fit <-
    function(object, parm, level=0.95,
             method=c("likelihood", "Wald"), ...)
{
    x <- object
    method <- match.arg(method)
    stopifnot(length(level) == 1L,
              is.numeric(level),
              level > 0,
              level < 1)
    ## test level.
    Data <- x$data
    if(method == "Wald") {
        if(is.null(x$vcov)) {
            ## Refit model to get vcov/std.err.:
            ctrl <- x$control
            ctrl$get.grad <- ctrl$get.vcov <- TRUE
            ctrl$test.args <- FALSE
            x <- dod_fit(same=Data[1, ], diff=Data[2, ], control=ctrl)
        }
        std.err <- NA_real_
        if(!is.null(x$vcov) && all(is.finite(x$vcov)))
            std.err <- sqrt(rev(diag(x$vcov))[1L])
        ci <- rep(NA_real_, 2L)
        if(is.finite(std.err)) {
            a <- 1 - level
            ci <- x$d.prime + qnorm(1 - a/2) * c(-1, 1) * std.err
            ci[1L] <- delimit(ci[1L], lower=0)
        }
    } else if(method == "likelihood") {
        ci <- get_ci(same=Data[1, ], diff=Data[2, ], alpha=1 - level,
                     dod.fit=x)
    } else
        stop("'method' not recognized")
    setNames(ci, nm=c("lower", "upper"))
}
