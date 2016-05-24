### <======================================================================>
"fit.ghypmv" <- function(data, lambda = 1, alpha.bar = 1,
                         mu = NULL, sigma = NULL, gamma = NULL,
                         opt.pars = c(lambda = T, alpha.bar = T, mu = T, sigma = T, gamma = !symmetric),
                         symmetric = F, standardize = F, nit = 2000, reltol = 1e-8, abstol = reltol * 10,
                         na.rm = F, silent = FALSE, save.data = T, trace = TRUE, ...)
{
    call <- match.call(expand.dots = TRUE)

    data <- .check.data(data = data, case = "mv", na.rm = na.rm)
    n <- nrow(data)
    d <- ncol(data)

    tmp.abar2chipsi <- .abar2chipsi(alpha.bar, lambda)
    chi <- tmp.abar2chipsi$chi
    psi <- tmp.abar2chipsi$psi

    ## check parameters of the mixing distribution for consistency
    .check.gig.pars(lambda, chi, psi)

    ## check opt.pars  for consistency
    opt.pars <- .check.opt.pars(opt.pars, symmetric)
    if(is.null(mu)){
        mu <- colMeans(data)
    }
    if(is.null(gamma) | symmetric){
        gamma <- rep(0, d)
    }
    if(symmetric){
        opt.pars["gamma"] <- FALSE
    }
    if(is.null(sigma)){
        sigma <- var(data)
    }

    ## check normal and skewness parameters  for consistency
    .check.norm.pars(mu, sigma, gamma, d)

    i.backup <- 1

    trace.pars <- list(alpha.bar = alpha.bar, lambda = lambda, mu = matrix(mu, ncol = d),
                       sigma = list(sigma), gamma = matrix(gamma, ncol = d))

    ## Internal function .fit.ghypmv allows reasonable error handling
    .fit.ghypmv <- function(data, lambda, alpha.bar, mu, sigma, gamma, opt.pars,
                            standardize, nit, reltol, abstol, silent, save.data, ...)
    {
        if(standardize){
            ## data will be standardized and initial values will be adapted
            tmp.mean <- colMeans(data)

            sigma.chol <- t(chol(solve(var(data))))

            data <- apply(data, MARGIN = 2, FUN = function(x){x - mean(x)})
            data <-  data %*% sigma.chol

            sigma <- t(sigma.chol) %*% sigma %*% sigma.chol
            gamma <- as.vector(sigma.chol %*% gamma)
            mu <- as.vector(sigma.chol %*% (mu - tmp.mean))
        }

        ## loglikelihood function of the inverse gamma distribution
        t.optfunc <- function(thepars, delta.sum, xi.sum, n.rows){
            nu <- -2 * .t.transform(thepars)
            term1 <- -n.rows * nu * log(nu/2 - 1)/2
            term2 <- (nu/2 + 1) * xi.sum + (nu/2 - 1) * delta.sum
            term3 <- n.rows * lgamma(nu/2)
            out <- term1 + term2 + term3
            return(out)
        }

        ## loglikelihood function of the gamma distribution
        vg.optfunc <- function(thepars, xi.sum, eta.sum, n.rows){
            thepars <- exp(thepars)
            term1 <- n.rows * (thepars * log(thepars) - lgamma(thepars))
            term2 <- (thepars - 1) * xi.sum - thepars * eta.sum
            out <- -(term1 + term2)
            return(out)
        }

        ## loglikelihood function of the generalized inverse gaussian distribution
        gig.optfunc <- function(thepars, mix.pars.fixed, delta.sum, eta.sum, xi.sum, n.rows)
        {
            out <- NA
            tmp.pars <- c(thepars, mix.pars.fixed)
            lambda <- tmp.pars["lambda"]
            alpha.bar <- exp(tmp.pars["alpha.bar"])
            tmp.abar2chipsi <- .abar2chipsi(alpha.bar, lambda)
            chi <- tmp.abar2chipsi$chi
            psi <- tmp.abar2chipsi$psi
            if (lambda < 0 & psi == 0){ # t
                out <-  t.optfunc(lambda, delta.sum, xi.sum, n.rows)
            }else if(lambda > 0 & chi == 0) { # VG
                out <-  vg.optfunc(lambda, xi.sum, eta.sum, n.rows)
            }else{                      # ghyp, hyp, NIG
                term1 <- (lambda - 1) * xi.sum
                term2 <- -chi * delta.sum/2
                term3 <- -psi * eta.sum/2
                term4 <- -n.rows * lambda * log(chi)/2 + n.rows * lambda * log(psi)/2 -
                    n.rows * .besselM3(lambda, sqrt(chi * psi), logvalue = TRUE)
                out <- -(term1 + term2 + term3 + term4)
            }
            return(out)
        }

        ## Initialize fitting loop
        i <- 0
        rel.closeness <- 100
        abs.closeness <- 100
        tmp.fit <- list(convergence = 0, message = NULL)

        tmp.abar2chipsi <- .abar2chipsi(alpha.bar, lambda)
        chi <- tmp.abar2chipsi$chi
        psi <- tmp.abar2chipsi$psi

        ll <-  sum(.dghypmv(data, lambda = lambda, chi = chi, psi = psi,
                            mu = mu, sigma = sigma, gamma = gamma, logvalue = TRUE))

        ## Start interations
        while ((abs.closeness > abstol) & (rel.closeness > reltol) & (i < nit)){
            i <- i + 1
            i.backup <<- i

            ##<------------------------ E-Step: EM update ------------------------------>
            ## The parameters mu, sigma and gamma become updated
            inv.sigma <- solve(sigma)
            Q <- mahalanobis(data, mu, inv.sigma, inverted = TRUE)
            Offset <- t(gamma) %*% inv.sigma %*% gamma
            delta <- Egig(lambda-d/2, Q+chi, psi+Offset, func = "1/x", check.pars = FALSE)
            delta.bar <- mean(delta)
            eta <- Egig(lambda-d/2, Q+chi, psi+Offset, func = "x", check.pars = FALSE)
            eta.bar <- mean(eta)
            delta.matrix <- matrix(delta, nrow = n, ncol = d, byrow = FALSE)
            if (opt.pars["gamma"]) {
                Xbar.matrix <- matrix(apply(data, 2, mean), nrow = n, ncol = d, byrow = TRUE)
                Xbar.matrix <- Xbar.matrix - data
                gamma <- apply(delta.matrix * Xbar.matrix, 2, sum) / (n * delta.bar * eta.bar - n)
            }
            if (opt.pars["mu"]) {
                mu <- (apply(delta.matrix * data, 2, sum)/n - gamma)/delta.bar
            }
            mu.matrix <- matrix(mu, nrow = n, ncol = d, byrow = TRUE)
            standardised <- data - mu.matrix
            tmp <- delta.matrix * standardised
            if (opt.pars["sigma"]) {
                sigma <- (t(tmp) %*% standardised)/n - outer(gamma, gamma) * eta.bar
            }

            ##<------------------------ M-Step: EM update ------------------------------>
            ## Maximise the conditional likelihood function and estimate lambda, chi, psi
            inv.sigma <- solve(sigma)
            Q <- mahalanobis(data, mu, inv.sigma, inverted = TRUE)
            Offset <- t(gamma) %*% inv.sigma %*% gamma

            xi.sum <- sum(Egig(lambda - d/2, Q + chi, psi + Offset, func = "logx", check.pars = FALSE))

            if(alpha.bar==0 & lambda > 0 & !opt.pars["alpha.bar"] & opt.pars["lambda"]){
                ##<------  VG case  ------>
                eta.sum <- sum(Egig(lambda - d/2, Q + chi, psi + Offset, func = "x", check.pars = FALSE))

                ## Supress warnings because if some parameters are not subject to optimization the
                ## problem can reduce to 1 dimension for which 'optim' is not optimal. Nonetheless,
                ## we use optim for the purpose to get a uniform list returned from one optimization routine.
                tmp.fit <- suppressWarnings(optim(log(lambda), vg.optfunc, eta.sum = eta.sum,
                                                  xi.sum = xi.sum, n.rows = n,...))
                lambda <- exp(tmp.fit$par)
            }else if(alpha.bar == 0 & lambda < 0 & !opt.pars["alpha.bar"] & opt.pars["lambda"]){
                ##<------  Student-t case  ------>
                delta.sum <- sum(Egig(lambda - d/2, Q + chi, psi + Offset, func = "1/x", check.pars = FALSE))

                tmp.fit <- suppressWarnings(optim(.inv.t.transform(lambda), t.optfunc, delta.sum = delta.sum,
                                                  xi.sum = xi.sum, n.rows = n,...))
                lambda <- .t.transform(tmp.fit$par)
            }else if(opt.pars["lambda"] | opt.pars["alpha.bar"]){
                ##<------  ghyp, hyp, NIG case  ------>
                delta.sum <- sum(Egig(lambda-d/2, Q+chi, psi+Offset, func="1/x", check.pars = FALSE))
                eta.sum <- sum(Egig(lambda-d/2, Q+chi, psi+Offset, func="x", check.pars = FALSE))

                mix.pars <- c(lambda = unname(lambda), alpha.bar = log(unname(alpha.bar)))
                opt.pars.mix <- opt.pars[c("lambda", "alpha.bar")]
                thepars <- mix.pars[opt.pars.mix]

                tmp.fit <- suppressWarnings(optim(thepars, gig.optfunc,
                                                  mix.pars.fixed = mix.pars[!opt.pars.mix],
                                                  delta.sum = delta.sum, eta.sum = eta.sum,
                                                  xi.sum = xi.sum, n.rows = n,...))

                lambda <- c(tmp.fit$par,mix.pars[!opt.pars.mix])["lambda"]
                alpha.bar <- exp(c(tmp.fit$par,mix.pars[!opt.pars.mix])["alpha.bar"])
            }
            tmp.abar2chipsi <- .abar2chipsi(alpha.bar, lambda)
            chi <- tmp.abar2chipsi$chi
            psi <- tmp.abar2chipsi$psi

            ## Test for convergence
            ll.old <- ll
            ll <-  sum(.dghypmv(data, lambda = lambda, chi = chi, psi = psi, mu = mu,
                                sigma = sigma, gamma = gamma, logvalue = TRUE))

            abs.closeness <- abs(ll - ll.old)
            rel.closeness <- abs((ll - ll.old)/ll.old)
            if(!is.finite(abs.closeness) | !is.finite(rel.closeness)){
                warning("fit.ghypmv: Loglikelihood is not finite! Iteration stoped!\n",
                        "Loglikelihood :", ll)
                break
            }

            ## Print result
            if(!silent){
                message <- paste("iter: ", i, "; rel.closeness: ", sprintf("% .6E", rel.closeness),
                                 "; log-likelihood: ", sprintf("% .6E", ll), "; alpha.bar: ", sprintf("% .6E", alpha.bar),
                                 "; lambda: ", sprintf("% .6E", lambda), sep = "")
                print(message)
            }
            ## Assign current values of optimization parameters
            ## to backups of the nesting environment.
            trace.pars$lambda <<- c(trace.pars$lambda, lambda)
            trace.pars$alpha.bar <<- c(trace.pars$alpha.bar, alpha.bar)
            trace.pars$mu <<- rbind(trace.pars$mu, mu)
            trace.pars$sigma <<- c(trace.pars$sigma, list(sigma))
            trace.pars$gamma <<- rbind(trace.pars$gamma, gamma)
        }
        ## END OF WHILE LOOP

        conv <- tmp.fit$convergence
        conv.type <- tmp.fit$message
        if(is.null(tmp.fit$message)){
            conv.type <- ""
        }else{
            conv.type <- paste("Message from 'optim':", tmp.fit$message)
        }

        converged <- FALSE
        if(i < nit & is.finite(rel.closeness) & is.finite(abs.closeness)){
            converged <- TRUE
        }

        if(standardize){
            inv.sigma.chol <- solve(sigma.chol)
            mu <- as.vector(inv.sigma.chol %*% mu + tmp.mean)
            sigma <- t(inv.sigma.chol) %*% sigma %*% inv.sigma.chol
            gamma <- as.vector(inv.sigma.chol %*% gamma)
            abar.chi.psi <- .abar2chipsi(alpha.bar = alpha.bar, lambda = lambda)
            ll <- sum(.dghypmv(x = data, lambda = lambda, chi = abar.chi.psi$chi,
                               psi = abar.chi.psi$psi, mu = mu, sigma = sigma,
                               gamma = gamma, logvalue = TRUE))
        }


        return(list(lambda = lambda, alpha.bar = alpha.bar, mu = mu, sigma = sigma, gamma = gamma,
                    llh = ll, n.iter = i, converged = converged, error.code = conv,
                    error.message = conv.type))

        ##  End of internal function .fit.ghypmv
    }



    save.fit <- try(.fit.ghypmv(data, lambda, alpha.bar, mu, sigma, gamma, opt.pars,
                                standardize, nit, reltol, abstol, silent, save.data, ...))


    if(class(save.fit) == "try-error"){
        lambda <- trace.pars$lambda[length(trace.pars$lambda)]
        alpha.bar <- trace.pars$alpha.bar[length(trace.pars$alpha.bar)]
        mu <- trace.pars$mu[nrow(trace.pars$mu), ]
        sigma <- trace.pars$sigma[[length(trace.pars$sigma)]]
        gamma <- trace.pars$gamma[nrow(trace.pars$gamma), ]
        llh <- as.numeric(NA)
        converged <- FALSE
        error.code <- 100
        error.message <- as.character(save.fit)
        n.iter <- i.backup
    }else{
        lambda <- save.fit$lambda
        alpha.bar <- save.fit$alpha.bar
        mu <- save.fit$mu
        sigma <- save.fit$sigma
        gamma <- save.fit$gamma
        llh <- save.fit$llh
        n.iter <- save.fit$n.iter
        converged <- save.fit$converged
        error.code <- save.fit$error.code
        error.message <- save.fit$error.message
    }

    if(!save.data){
        data <- NULL
    }

    if(!trace){
        trace.pars <- list()
    }

    nbr.fitted.params <- unname(sum(opt.pars[c("alpha.bar","lambda")]) +
                                d * sum(opt.pars[c("mu","gamma")]) +
                                d/2 * (d + 1) * opt.pars[c("sigma")])
    aic <- -2 * llh + 2 * nbr.fitted.params

    ghyp.object <- ghyp(lambda = lambda, alpha.bar = alpha.bar,
                        mu = mu, sigma = sigma, gamma = gamma, data = data)

    ghyp.object@parametrization <- "alpha.bar"
    ghyp.object@call <- call

    return(.fit.ghyp(ghyp.object, llh = llh, n.iter = n.iter, converged = converged,
                     error.code = error.code, error.message = error.message,
                     fitted.params = opt.pars, aic = aic, trace.pars = trace.pars))
}
### <---------------------------------------------------------------------->


### <======================================================================>
"fit.ghypuv" <- function(data, lambda = 1, alpha.bar = 0.5, mu = median(data),
                         sigma = mad(data), gamma = 0,
                         opt.pars = c(lambda = T, alpha.bar = T, mu = T, sigma = T, gamma = !symmetric),
                         symmetric = F, standardize = F, save.data = T, na.rm = T, silent = FALSE, ...)
{
    call <- match.call(expand.dots = TRUE)

    ## Check input
    opt.pars <- .check.opt.pars(opt.pars, symmetric)
    data <- .check.data(data = data, case = "uv", na.rm = na.rm)
    if(standardize){
        ## data will be standardized and initial values will be adapted
        tmp.mean <- mean(data)
        tmp.sd <- sd(data)
        data.backup <- data
        data <- (data - tmp.mean)/tmp.sd
        mu <- (mu - tmp.mean)/tmp.sd
        sigma <- sigma/tmp.sd
        gamma <- gamma/tmp.sd
    }

    .check.norm.pars(mu, sigma, gamma, 1)

    tmp.abar2chipsi <- .abar2chipsi(alpha.bar, lambda)
    chi <- tmp.abar2chipsi$chi
    psi <- tmp.abar2chipsi$psi

    ## check parameters of the mixing distribution for consistency
    .check.gig.pars(lambda, chi, psi)

    var.names <- c("lambda", "alpha.bar", "mu", "sigma", "gamma")

    vars <- unname(c(lambda, alpha.bar, mu, sigma, gamma))
    names(vars) <- var.names

    if(alpha.bar == 0 & !opt.pars["alpha.bar"]){
        if(lambda > 0){                 # VG case
            transform <- c("exp", "exp")
            inv.transform <- c("log", "log")
        }else{                          # Student-t case
            transform <- c(".t.transform", "exp")
            inv.transform <- c(".inv.t.transform", "log")
        }
        names(transform) <- names(inv.transform) <- c("lambda", "sigma")
    }else{                              # ghyp, NIG, hyp case
        transform <- c("exp", "exp")
        inv.transform <- c("log", "log")
        names(transform) <- names(inv.transform) <- c("alpha.bar", "sigma")
    }



    ## Inverse transformation of the initial parameter values
    for(nam in intersect(names(opt.pars[opt.pars]), names(transform))) {
        vars[nam] <- do.call(inv.transform[nam], list(vars[nam]))
    }
    tmp.fit <- .mle.default(data = data, pdf = ".dghypuv",
                            vars = vars, opt.pars = opt.pars,
                            transform = transform, se = T,
                            silent = silent, ...)

    if(tmp.fit$convergence != 0){
        converged <- FALSE
    }else{
        converged <- TRUE
    }

    if(is.null(tmp.fit$message)){
        conv.type <- ""
    }else{
        conv.type <- paste("Message from 'optim':", tmp.fit$message)
    }

    if(standardize){
        data <- data.backup
        tmp.fit$par.ests["mu"] <- tmp.fit$par.ests["mu"] * tmp.sd + tmp.mean
        tmp.fit$par.ests["sigma"] <- tmp.fit$par.ests["sigma"] * tmp.sd
        tmp.fit$par.ests["gamma"] <- tmp.fit$par.ests["gamma"] * tmp.sd
        tmp.llh <- try(sum(.dghypuv(x = data, lambda = tmp.fit$par.ests["lambda"],
                                    alpha.bar = tmp.fit$par.ests["alpha.bar"],
                                    mu = tmp.fit$par.ests["mu"],
                                    sigma = tmp.fit$par.ests["sigma"],
                                    gamma = tmp.fit$par.ests["gamma"],
                                    logvalue = TRUE)))
        if(class(tmp.llh) == "try-error"){
            warning("Error occured during renormalization! Log-likelihood set to zero!\n")
            tmp.fit$ll.max <- as.numeric(NA)
        }else{
            tmp.fit$ll.max <- tmp.llh
        }
        tmp.fit$parameter.variance <- matrix(NA)
    }

    nbr.fitted.params <- unname(sum(opt.pars))
    aic <- -2 * tmp.fit$ll.max + 2 * nbr.fitted.params

    if(!save.data){
        data <- NULL
    }

    ghyp.object <- ghyp(lambda = tmp.fit$par.ests["lambda"],
                        alpha.bar = tmp.fit$par.ests["alpha.bar"],
                        mu = tmp.fit$par.ests["mu"],
                        sigma = tmp.fit$par.ests["sigma"],
                        gamma = tmp.fit$par.ests["gamma"],
                        data = data)

    ghyp.object@call <- call
    return(.fit.ghyp(ghyp.object, llh = tmp.fit$ll.max, n.iter = tmp.fit$n.iter,
                    converged = converged,
                    error.code = tmp.fit$convergence, error.message = conv.type,
                    parameter.variance = tmp.fit$parameter.variance,
                    fitted.params = opt.pars, aic = aic, trace.pars = as.list(as.data.frame(tmp.fit$trace.pars))))
}
### <---------------------------------------------------------------------->


### <======================================================================>
"stepAIC.ghyp" <- function(data, dist = c("ghyp", "hyp", "NIG", "VG", "t", "gauss"),
                           symmetric = NULL, ...)
{
    call <- match.call(expand.dots = TRUE)
    dist <- match.arg(dist, several.ok = TRUE)

    with.gauss <- "gauss" %in% dist

    dist <- dist[!(dist == "gauss")]

    type <- "uv"
    tmp.data <- try(.check.data(data, case = "uv", na.rm = T, fit = TRUE), silent = TRUE)
    if(class(tmp.data) == "try-error"){
        tmp.data <- try(.check.data(data, case = "mv", na.rm = T, fit = TRUE), silent = TRUE)
        type <- "mv"
        if(class(tmp.data) == "try-error"){
            stop("Invalid data!")
        }
    }
    if(is.null(symmetric)){
        symm <- c(FALSE, TRUE)
    }else{
        if(symmetric == FALSE){
            symm <- FALSE
        }else if(symmetric == TRUE){
            symm <- TRUE
        }else{
            stop("Argument 'asymmetric' must be either 'NULL', 'TRUE' or 'FALSE'!")
        }
    }
    nbr.fits <- length(dist) * length(symm)
    function.names <- paste("fit.", dist, type, sep = "")
    fitted.objects <- vector(mode = "list", nbr.fits)

    fit.info <- data.frame(model = rep(dist, length(symm)),
                           symmetric = rep(symm, each = length(dist)),
                           lambda = rep(NA, nbr.fits),
                           alpha.bar = rep(NA, nbr.fits),
                           aic = rep(NA, nbr.fits),
                           llh = rep(NA, nbr.fits),
                           converged = rep(NA, nbr.fits),
                           n.iter = rep(NA, nbr.fits),
                           stringsAsFactors = FALSE)

    ## In the univariate case return a data.frame with each ghyp parameter
    if(type=="uv"){
        uv.params <- data.frame(mu = rep(NA, nbr.fits),
                                sigma = rep(NA, nbr.fits),
                                gamma = rep(NA, nbr.fits),
                                stringsAsFactors = FALSE)
    }

    for(j in 1:length(symm)){
        for(i in 1:length(function.names)){
            if(symm[j]){
                cat("Currently fitting: symmetric", dist[i], "\n")
            }else{
                cat("Currently fitting: asymmetric", dist[i], "\n")
            }
            call.args <- c(list(data = tmp.data), list(...), list(symmetric = symm[j]))

            tmp.fit <- do.call(function.names[i], call.args)

            tmp.fit@call <- call

            tmp.params <- coef(tmp.fit, type = "alpha.bar")
            tmp.fit.info <- ghyp.fit.info(tmp.fit)

            if(type=="uv"){
                uv.params[(j - 1) * length(dist) + i, ] <- data.frame(tmp.params$mu,
                                                                      tmp.params$sigma,
                                                                      tmp.params$gamma,
                                                                      stringsAsFactors = FALSE)
            }

            tmp.result <- data.frame(dist[i], symm[j],
                                     tmp.params$lambda, tmp.params$alpha.bar,
                                     tmp.fit.info$aic, tmp.fit.info$logLikelihood,
                                     tmp.fit.info$converged, tmp.fit.info$n.iter,
                                     stringsAsFactors = FALSE)

            fit.info[(j - 1) * length(dist) + i, ] <- tmp.result
            fitted.objects[[(j - 1) * length(dist) + i]] <- tmp.fit
        }
    }

    if(with.gauss){
        if(!is.null(list(...)$save.data)){
            save.data <- list(...)$save.data
        }else{
            save.data <- T
        }
        if(!is.null(list(...)$na.rm)){
            na.rm <- list(...)$na.rm
        }else{
            na.rm <- T
        }
        cat("Currently fitting: gauss\n")
        if(type == "uv"){
            gauss.fit <- fit.gaussuv(tmp.data, save.data = save.data, na.rm = na.rm)
            gauss.uv.params <- data.frame(gauss.fit@mu, gauss.fit@sigma, 0, stringsAsFactors = FALSE)
            colnames(gauss.uv.params) <- colnames(uv.params)
            uv.params <- rbind(uv.params, gauss.uv.params)
        }else{
            gauss.fit <- fit.gaussmv(tmp.data, save.data = save.data, na.rm = na.rm)
        }

        gauss.fit.info <- data.frame("gauss", TRUE, NA, Inf, AIC(gauss.fit),
                                     logLik(gauss.fit), TRUE, 0, stringsAsFactors = FALSE)

        colnames(gauss.fit.info) <- colnames(fit.info)
        fit.info <- rbind(fit.info, gauss.fit.info)

        fitted.objects <- c(fitted.objects, gauss.fit)
    }
    if(type == "uv"){
        fit.info <- cbind(fit.info[, 1:4], uv.params, fit.info[, 5:8])
    }

    idx <- which(fit.info$aic == min(fit.info$aic, na.rm = TRUE))
    if(length(idx) > 1){
        warning("Several AIC minima observed; the first minima will be returned!")
        idx <- idx[1]
    }
    best.model <- fitted.objects[[idx]]
    fit.info <- fit.info[order(fit.info$aic, na.last = TRUE), ]
    return(list(best.model = best.model, all.models = fitted.objects, fit.table = fit.info))
}
### <---------------------------------------------------------------------->
