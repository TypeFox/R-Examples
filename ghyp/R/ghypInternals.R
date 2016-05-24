### <======================================================================>
".abar2chipsi" <- function(alpha.bar, lambda, eps = .Machine$double.eps)
{
    if(alpha.bar < 0){
        stop("alpha.bar must be non-negative.")
    }
    if(alpha.bar > eps){
        if(lambda >= 0){
            psi = alpha.bar * besselK(alpha.bar, lambda + 1, expon.scaled = TRUE)/
                besselK(alpha.bar, lambda, expon.scaled = TRUE)
            if(is.na(psi)){
                psi <- 200
                ## warning( paste( "Overflow in besselK, alpha.bar = ",
                ##               alpha.bar, "; lambda = ", lambda) )
            }
            chi <- alpha.bar^2 / psi
        }else{
            chi <- alpha.bar * besselK(alpha.bar, lambda, expon.scaled = TRUE)/
                besselK(alpha.bar, lambda + 1, expon.scaled = TRUE)
            if(is.na(chi)){
                chi <- 200
                ## warning( paste( "Overflow in besselK, alpha.bar = ",
                ##               alpha.bar, "; lambda = ", lambda) )
            }
            psi <- alpha.bar^2 / chi
        }
    }else{
        if(lambda > 0){                 # VG
            chi <- 0
            psi <- 2 * lambda
        }else if(lambda < 0){           # Student-t
            psi <- 0
            chi <- -2 * (lambda + 1)
        }else{
            stop("Forbidden combination of parameter values");
        }
    }
    list("lambda" = unname(lambda), "chi" = unname(chi), "psi" = unname(psi))
}
### <---------------------------------------------------------------------->



### <======================================================================>
".besselM3" <- function(lambda = 9/2, x = 2, logvalue = FALSE)
{
    if(all(abs(lambda) == 0.5)){
        ## Simplified expression in case lambda == 0.5
        if (!logvalue){
            res <- sqrt(pi/(2 * x)) * exp(-x)
        }else{
            res <- 0.5 * log(pi/(2 * x)) - x
        }
    }else{
        if (!logvalue){
            res <- besselK(x, lambda)
        }else{
            res <- log( besselK(x, lambda, expon.scaled = TRUE) ) - x
        }
    }
    return(res)
}
### <---------------------------------------------------------------------->



### <======================================================================>
".check.data" <- function(data, case = c("uv", "mv"), na.rm = TRUE,
                         fit = TRUE, dim = NULL)
{
    case <- match.arg(case)

    if(case == "mv"){                   # MULTIVARIATE
        data <- as.matrix(data)
        ##    if(ncol(data) < 2){
        ##      stop("'data' must have more than one column. ",
        ##           "Try the univariate case!")
        ##    }
        if(any(dim(data) == 1)){
            data <- matrix(data, nrow = 1)
        }
        na.idx <- apply(is.na(data), 1, any)
        if(na.rm){
            if(any(na.idx)){
                if(all(na.idx)){
                    stop("Sample contains only columns with NA's!")
                }else{
                    warning(sum(na.idx)," NA observations removed")
                }
                data <- data[!na.idx, ]
            }
        }else{
            if(all(na.idx)){
                stop("Sample contains only columns with NA's!")
            }
        }
        if(nrow(data) < ncol(data) & fit){
            stop("'data' must not have more columns than rows.")
        }
    }else{
        if(!is.vector(data)){           # UNIVARIATE
            data <- as.matrix(data)
            if(ncol(data) > 1){
                stop("'data' must have only one column!")
            }
            data <- as.vector(data)
        }

        na.idx <- is.na(data)
        if(na.rm){
            if(any(na.idx)){
                if(all(na.idx)){
                    stop("Sample contains only NA's!")
                }else{
                    warning(sum(na.idx), " NA observations removed")
                }
            }
            data <- data[!na.idx]
        }else{
            if(all(na.idx)){
                stop("Sample contains only NA's!")
            }
        }

        if(fit & length(data) < 5){
            stop("Length of 'data' must be at least 5.")
        }
    }
    if(!is.numeric(data)){
        stop("'data' must be numeric!")
    }
    return(data)
}
### <---------------------------------------------------------------------->



### <======================================================================>
".check.gig.pars" <- function(lambda, chi, psi)
{
    ## Vectorized check
    "internal.check.gig.pars" <- function(x){
        lambda <- x[1]
        chi <- x[2]
        psi <- x[3]
        if(lambda < 0 & (chi <= 0 | psi < 0)){
            stop("If lambda < 0: chi must be > 0 and psi must be >= 0! \n",
                 "lambda = ", lambda, ";   chi = ", chi, ";   psi = ", psi, "\n")
        }
        if(lambda == 0 & (chi <= 0 | psi <= 0)){
            stop("If lambda == 0: chi must be > 0 and psi must be > 0! \n",
                 "lambda = ", lambda, ";   chi = ", chi, ";   psi = ", psi, "\n")
        }
        if(lambda > 0 & (chi < 0 | psi <= 0)){
            stop("If lambda > 0: chi must be >= 0 and psi must be > 0! \n",
                 "lambda = ", lambda, ";   chi = ", chi, ";   psi = ", psi, "\n")
        }
    }

    params <- suppressWarnings(cbind(lambda, chi, psi))
    apply(params, 1, internal.check.gig.pars)
}
### <---------------------------------------------------------------------->



### <======================================================================>
".check.norm.pars" <- function(mu, sigma, gamma, dimension)
{
    if(length(mu) != dimension){
        stop("Parameter 'mu' must be of length ", dimension, "!")
    }
    if(length(gamma) != dimension){
        stop("Parameter 'gamma' must be of length ", dimension, "!")
    }
    if(dimension > 1){                  # MULTIVARIATE
        if(!is.matrix(sigma)){
            stop("'sigma' must be a quadratic matrix with dimension ",
                 dimension, " x ", dimension, "!")
        }
        if(nrow(sigma) != dimension | ncol(sigma) != dimension){
            stop("Matrix 'sigma' must be quadratic with dimension ",
                 dimension, " x ", dimension, "!")
        }
    }else{                              # UNIVARIATE
        if(length(sigma) != dimension){
            stop("Parameter 'sigma' must be a scalar!")
        }
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
".check.opt.pars" <- function(opt.pars, symmetric)
{
    default.names <- c("lambda", "alpha.bar", "mu", "sigma", "gamma")
    ## There are 3 possibilities:
    ## (1) opt.pars are not named: They must be of length 5
    ##     and are interpreted in the same order as 'default.names'
    ## (2) opt.pars are named and are in 'default.names':
    ##     (*) Any unknown names are dropped
    ##     (*) opt.pars is reordered
    ## (3) If opt.pars with unknown names are passed an error occurs
    if(is.null(names(opt.pars))){
        if(length(opt.pars) == 5){
            new.opt.pars <- c(lambda = opt.pars[1], alpha.bar = opt.pars[2],
                              mu = opt.pars[3], sigma = opt.pars[4],
                              gamma = opt.pars[5])

            if(new.opt.pars["gamma"] & symmetric){
                warning("Clash: symmetric == TRUE while opt.pars[5] == TRUE!\n",
                        "A symmetric model will be fitted (i.e. opt.pars[5] <- FALSE)!\n")
            }
        }else{
            stop("If 'opt.pars' is not named it must have length 5!\n",
                 "The order is lambda, alpha.bar, mu, sigma, gamma.\n")
        }
    }else{
        ## In case a named argument 'symmetric' was submitted and
        ## 'opt.pars' was not submited, the name corresponding to the
        ## 'gamma' parameter in 'opt.pars' is of the structure
        ## 'gamma.xx'. Here, set it back to 'gamma'.
        if(length(grep("gamma", names(opt.pars))) > 0){
            names(opt.pars)[grep("gamma", names(opt.pars))] <- "gamma"
        }

        if(all(default.names %in% names(opt.pars))){
            if(length(opt.pars) != 5){
                warning("The following names were dropped:\n",
                        paste(names(opt.pars)[!(names(opt.pars) %in% default.names)],
                              collapse = ", "))
            }
            new.opt.pars <- c(lambda = unname(opt.pars["lambda"]),
                              alpha.bar = unname(opt.pars["alpha.bar"]),
                              mu = unname(opt.pars["mu"]),
                              sigma = unname(opt.pars["sigma"]),
                              gamma = unname(opt.pars["gamma"]))
            if(new.opt.pars["gamma"] & symmetric){
                warning("Clash: symmetric == TRUE while opt.pars['gamma'] == TRUE!\n",
                        "A symmetric model will be fitted (i.e. opt.pars['gamma'] <- FALSE)!\n")
            }
        }else if(!any(default.names %in% names(opt.pars))){
            stop("The names '", paste(names(opt.pars), collapse = "', '"),
                 "' do not match the required names lambda, alpha.bar, mu, sigma, gamma.\n")
        }else if(!all(names(opt.pars) %in% default.names)){
            stop("The names '", paste(names(opt.pars)[!(names(opt.pars) %in% default.names)],
                                      collapse = "', '"),
                 "' do not match the required names lambda, alpha.bar, mu, sigma, gamma.\n")
        }else{
            new.opt.pars <- logical(0)
            if("lambda" %in% names(opt.pars)){
                new.opt.pars <-  c(new.opt.pars, lambda = unname(opt.pars["lambda"]))
            }else{
                new.opt.pars <-  c(new.opt.pars, lambda = TRUE)
            }
            if("alpha.bar" %in% names(opt.pars)){
                new.opt.pars <-  c(new.opt.pars, alpha.bar = unname(opt.pars["alpha.bar"]))
            }else{
                new.opt.pars <-  c(new.opt.pars, alpha.bar = TRUE)
            }
            if("mu" %in% names(opt.pars)){
                new.opt.pars <-  c(new.opt.pars, mu = unname(opt.pars["mu"]))
            }else{
                new.opt.pars <-  c(new.opt.pars, mu = TRUE)
            }
            if("sigma" %in% names(opt.pars)){
                new.opt.pars <-  c(new.opt.pars, sigma = unname(opt.pars["sigma"]))
            }else{
                new.opt.pars <-  c(new.opt.pars, sigma = TRUE)
            }
            if("gamma" %in% names(opt.pars)){
                new.opt.pars <-  c(new.opt.pars, gamma = unname(opt.pars["gamma"]))
                if(new.opt.pars["gamma"] & symmetric){
                warning("Clash: symmetric == TRUE while opt.pars['gamma'] == TRUE!\n",
                        "A symmetric model will be fitted (i.e. opt.pars['gamma'] <- FALSE)!\n")
                }
            }else{
                new.opt.pars <-  c(new.opt.pars, gamma = TRUE)
            }
        }
    }
    if(symmetric){
        new.opt.pars["gamma"] <- FALSE
    }
    return(new.opt.pars)
}
### <---------------------------------------------------------------------->



### <======================================================================>
".fit.ghyp" <- function(object, llh = 0, n.iter = 0, converged = FALSE, error.code = 0,
                       error.message = "", parameter.variance, fitted.params, aic,
                       trace.pars = list())
{
    if(missing(parameter.variance)){
        parameter.variance <- matrix(0)
    }

    return(new("mle.ghyp", call = object@call, lambda = object@lambda, chi = object@chi,
               psi = object@psi, alpha.bar = object@alpha.bar, mu = object@mu,
               sigma = object@sigma, gamma = object@gamma, model = object@model,
               dimension = object@dimension, expected.value = object@expected.value,
               variance = object@variance, data = object@data,
               parametrization = object@parametrization, llh =llh,
               n.iter = n.iter, converged = converged, error.code = error.code,
               error.message = error.message,
               parameter.variance = parameter.variance,
               fitted.params = fitted.params, aic = aic,
               trace.pars = trace.pars))
}
### <---------------------------------------------------------------------->



### <======================================================================>
".get.stepAIC.ghyp" <- function(stepAIC.obj,
                               dist = c("ghyp", "hyp", "NIG", "VG", "t", "gauss"),
                               symmetric = FALSE)
{
    dist <- match.arg(dist, several.ok = FALSE)

    if(is.null(symmetric) || is.na(symmetric)){
        stop("'symmetric' must be either 'TRUE' or 'FALSE'!")
    }

    if(dist == "gauss" && !symmetric){
        stop("Asymmetric 'gauss' does not exist!")
    }

    table.idx <- which(stepAIC.obj$fit.table$model == dist &
                       stepAIC.obj$fit.table$symmetric == symmetric)

    if(length(table.idx) == 0){
        stop("Model '", dist, "' for symmetric = ", symmetric, " not found!")
    }

    model.idx <- as.numeric(rownames(stepAIC.obj$fit.table)[table.idx])

    return(list(model = stepAIC.obj$all.models[[model.idx]],
                fit.info = stepAIC.obj$fit.table[table.idx, ]))
}
### <---------------------------------------------------------------------->



### <======================================================================>
".ghyp.model" <- function(lambda, chi, psi, gamma)
{
    d <- length(gamma)
    if(chi == Inf && psi == Inf){
        ghyp.case.long.skew <- "Gaussian"
        ghyp.case.long <- "Gaussian"
        ghyp.case.short.skew <- "Gauss"
        ghyp.case.short <- "Gauss"
    }else{
        if(lambda == (d + 1)/2 & chi > 0 & psi > 0){
            ghyp.case.long <- "Hyperbolic"
            ghyp.case.short <- "hyp"
        }else if(lambda == -0.5 & chi > 0 & psi > 0){
            ghyp.case.long <- "Normal Inverse Gaussian"
            ghyp.case.short <- "NIG"
        }else if(lambda > 0 & chi == 0){
            ghyp.case.long <- "Variance Gamma"
            ghyp.case.short <- "VG"
        }else if(lambda < 0 & psi == 0){
            ghyp.case.long <- "Student-t"
            ghyp.case.short <- "t"
        }else{
            ghyp.case.long <- "Generalized Hyperbolic"
            ghyp.case.short <- "ghyp"
        }
        if(sum(abs(gamma)) == 0){
            ghyp.case.long.skew <- paste("Symmetric", ghyp.case.long, sep = " ")
            ghyp.case.short.skew <- paste("Symm", ghyp.case.short, sep = " ")
        }else{
            ghyp.case.long.skew <- paste("Asymmetric", ghyp.case.long, sep = " ")
            ghyp.case.short.skew <- paste("Asymm", ghyp.case.short, sep = " ")
        }
    }
    return(unname(c(ghyp.case.long.skew, ghyp.case.long,
                    ghyp.case.short.skew, ghyp.case.short)))
}
### <---------------------------------------------------------------------->



### <======================================================================>
".integrate.moment.ghypuv" <- function(x, moment = 1, ...)
{
    .dghypuv(x, ...) * x^moment
}
### <---------------------------------------------------------------------->



### <======================================================================>
".integrate.moment.gig" <- function(x,moment=1,...)
{
    dgig(x,...) * x^moment
}
### <---------------------------------------------------------------------->



### <======================================================================>
".dghypuv" <- function(x, lambda = 1, chi = 1, psi = 1, alpha.bar = NULL,
                             mu = 1, sigma = 1, gamma = 0, logvalue = FALSE)
{
    sigma <- as.vector(sigma)
    if(!is.null(alpha.bar)){
        tmp.abar2chipsi <- .abar2chipsi(alpha.bar, lambda)
        chi <- tmp.abar2chipsi$chi
        psi <- tmp.abar2chipsi$psi
    }
    Q <- ((x - mu)/sigma)^2
    d <- 1
    n <- length(x)
    inv.sigma <- 1/sigma^2

    if(gamma == 0){
        symm <- TRUE
        skewness.scaled <- 0
        skewness.norm <- 0
    } else {
        symm <- FALSE
        skewness.scaled <- (x - mu) * (inv.sigma * gamma)
        skewness.norm <- gamma^2 * inv.sigma
    }
    out <- NA
    if (psi == 0){
        lambda.min.0.5 <- lambda - 0.5
        if(symm){                       # Symmetric Student-t
##             nu <- -2 * lambda
##             sigma.t <- sqrt((nu - 2) / nu) * sigma
##             out <- dt((x - mu) / sigma.t, df = nu, log = TRUE) - log(sigma.t)

            interm <- chi + Q

            log.const.top <- -lambda * log(chi) + lgamma(-lambda.min.0.5)
            log.const.bottom <- 0.5 * log(pi) + log(sigma) + lgamma(-lambda)
            log.top <- lambda.min.0.5 * log(interm)

            out <- log.const.top + log.top - log.const.bottom
        }else{                          # Asymmetric Student-t
            interm <- sqrt((chi + Q) * skewness.norm)

            log.const.top <- -lambda * log(chi) - lambda.min.0.5 * log(skewness.norm)
            log.const.bottom <- 0.5 * log(2 * pi) + log(sigma) + lgamma(-lambda) - (lambda + 1) * log(2)

            log.top <- .besselM3(lambda.min.0.5, interm, logvalue = TRUE) + skewness.scaled
            log.bottom <- -lambda.min.0.5 * log(interm)

            out <- log.const.top + log.top - log.const.bottom - log.bottom
        }
    }else if(psi > 0){
        if(chi > 0){                    # ghyp, hyp and NIG (symmetric and asymmetric)
            log.top <- .besselM3((lambda - 0.5), sqrt((psi + skewness.norm) * (chi + Q)),logvalue = TRUE) + skewness.scaled
            log.bottom <- (0.5 - lambda) * log(sqrt((psi + skewness.norm) * (chi + Q)))
            log.const.top <- -lambda/2 * log(psi * chi) + 0.5 * log(psi) +
                (0.5 - lambda ) * log(1 + skewness.norm/psi)
            log.const.bottom <- 0.5 * log(2 * pi) +
                .besselM3(lambda, sqrt(chi * psi),logvalue = TRUE) + log(sigma)
            out <- log.const.top + log.top - log.const.bottom - log.bottom
        }else if(chi == 0){             # Variance gamma (symmetric and asymmetric)
            ## Density contains singularities for Q -> 0. Three cases depending on
            ## lambda are catched:

            ## any(Q < eps) & lambda >  0.5: Interpolate with splines
            ## any(Q < eps) & lambda <  0.5: Set Q[Q < eps] to eps
            ## any(Q < eps) & lambda == 0.5: Set Q[Q < eps] to NA

            ##<---------------  Internal function vg.density  ------------------>
            vg.density <- function(Q, skewness.norm, skewness.scaled, lambda, psi, sigma)
            {
                log.top <- .besselM3((lambda - 0.5), sqrt((psi + skewness.norm) * Q),logvalue = TRUE) + skewness.scaled
                log.bottom <- (0.5 - lambda) * log(sqrt((psi + skewness.norm) * Q))

                log.const.top <- log(psi) * lambda + (1 - lambda) * log(2) +
                    (0.5 - lambda) * log(psi + skewness.norm)

                log.const.bottom <- 0.5 * log(2 * pi) + lgamma(lambda) + log(sigma)
                return(log.const.top + log.top - log.const.bottom - log.bottom)
            }
            ##<-------------  End of internal function vg.density  ------------>

            eps <- .Machine$double.eps
            ## Observations that are close to mu were kept at a minimum magnitude
            if(any(Q < eps)){
                ## If lambda == 0.5 * dimension, there is another singularity.
                if(lambda == 0.5){
                    warning("NA's generated in .dghypuv: ",
                            "Some standardized observations are close to 0 ",
                            "and lambda is close to 0.5!")
                    if(gamma == 0){
                        tmp.skewness.scaled <- 0
                    }else{
                        tmp.skewness.scaled <- skewness.scaled[Q >= eps]
                    }
                    out <- rep(NA, length(Q))
                    out[Q >= eps] <- vg.density(Q[Q >= eps], skewness.norm, tmp.skewness.scaled, lambda, psi, sigma)
                }else if(lambda > 0.5){
                    message("Singularity (x-mu)==0: Interpolate with splines.")

                    ##<----------  Internal function vg.density.singular  -------------->
                    vg.density.singular <- function(skewness.norm, lambda, psi, sigma)
                    {
                        log.const.top <- log(psi) * lambda + (0.5 - lambda) * log(psi + skewness.norm) +
                            lgamma(lambda - 0.5)

                        log.const.bottom <- log(2) + 0.5 * log(pi) + lgamma(lambda) + log(sigma)
                        return(log.const.top - log.const.bottom)
                    }
                    ##<------  End of internal function vg.density.singular  ----------->

                    out <- rep(0, length(Q))

                    if(gamma == 0){
                        tmp.skewness.scaled <- 0
                    }else{              # Compute observations > eps as usual
                        tmp.skewness.scaled <- skewness.scaled[Q >= eps]
                    }

                    out[Q >= eps] <- vg.density(Q[Q >= eps], skewness.norm, tmp.skewness.scaled, lambda, psi, sigma)

                    ## Interpolate all observations < eps
                    ## x points
                    tmp.x <- c(-2, -1, 1, 2) * sqrt(eps) * sigma + mu
                    spline.x <- c(tmp.x[1:2], mu, tmp.x[3:4])

                    ## y points
                    if(gamma == 0){
                        tmp.skewness.scaled <- 0
                    }else{
                        tmp.skewness.scaled <- (tmp.x - mu) * inv.sigma * gamma
                    }

                    tmp.Q <- ((tmp.x - mu)/sigma)^2
                    tmp.density <-  vg.density(tmp.Q, skewness.norm, tmp.skewness.scaled, lambda, psi, sigma)

                    spline.y <- c(tmp.density[1:2], vg.density.singular(skewness.norm, lambda, psi, sigma),
                                  tmp.density[3:4])

                    vg.density.interp <- splinefun(spline.x, spline.y)
                    out[Q < eps] <- vg.density.interp(x[Q < eps])
                }else{
                    Q[Q < eps] <- eps
                    warning("Singularity: Some standardized observations are close to 0 (< ",
                            sprintf("% .6E", eps), ")! Observations set to ",
                            sprintf("% .6E", eps), ".", immediate. = TRUE)

                    out <- vg.density(Q, skewness.norm, skewness.scaled, lambda, psi, sigma)
                }
            }else{
                out <- vg.density(Q, skewness.norm, skewness.scaled, lambda, psi, sigma)
            }
        }
        else out <- NA
    }
    if (!logvalue){
        out <- exp(out)
    }
    return(out)
}
### <---------------------------------------------------------------------->



### <======================================================================>
".dghypmv" <- function(x, lambda, chi, psi, mu, sigma, gamma, logvalue = FALSE)
{
    ## Density of a multivariate generalized hyperbolic distribution.
    ## Covers all special cases as well.
    d <- dim(x)[2]
    n <- dim(x)[1]
    det.sigma <- det(sigma)
    inv.sigma <- solve(sigma)
    Q <- mahalanobis(x, mu, inv.sigma, inverted = TRUE)

    if (sum(abs(gamma)) == 0){
        symm <- TRUE
        skewness.scaled <- 0
        skewness.norm <- 0
    } else {
        symm <- FALSE
        skewness.scaled <- as.vector((as.matrix(x) - matrix(mu, nrow = n, ncol = d, byrow = T)) %*%
                                     (inv.sigma %*% gamma))
        skewness.norm <- t(gamma) %*% inv.sigma %*% gamma
    }
    out <- NA
    if (psi == 0){
        lambda.min.d.2 <- lambda - d / 2
        if(symm){                       # Symmetric Student-t
            interm <- chi + Q

            log.const.top <- -lambda * log(chi) + lgamma(-lambda.min.d.2)
            log.const.bottom <- d / 2 * log(pi) + 0.5 * log(det.sigma) + lgamma(-lambda)
            log.top <- lambda.min.d.2 * log(interm)

            out <- log.const.top + log.top - log.const.bottom

        }else{                          # Asymmetric Student-t
            interm <- sqrt((chi + Q) * skewness.norm)

            log.const.top <- -lambda * log(chi) - lambda.min.d.2 * log(skewness.norm)
            log.const.bottom <- d / 2 * log(2 * pi) + 0.5 * log(det.sigma) +
                lgamma(-lambda) - (lambda + 1) * log(2)
            log.top <- .besselM3(lambda.min.d.2, interm, logvalue = TRUE) + skewness.scaled
            log.bottom <- -lambda.min.d.2 * log(interm)

            out <- log.const.top + log.top - log.const.bottom - log.bottom

        }
    }
    else if (psi > 0){
        if (chi > 0){                   # ghyp, hyp and NIG (symmetric and asymmetric)
            log.top <- .besselM3((lambda - d/2), sqrt((psi + skewness.norm) * (chi + Q)),
                                logvalue = T) + skewness.scaled
            log.bottom <- (d/2 - lambda) * log(sqrt((psi + skewness.norm) * (chi + Q)))
            log.const.top <- -lambda/2 * log(psi * chi) + (d/2) * log(psi) +
                (d/2 - lambda ) * log(1 + skewness.norm/psi)
            log.const.bottom <- d/2 * log(2 * pi) +
                .besselM3(lambda, sqrt(chi * psi), logvalue = T) + 0.5 * log(det.sigma)
            out <- log.const.top + log.top - log.const.bottom - log.bottom
        }
        else if (chi == 0){             # Variance gamma (symmetric and asymmetric)
            eps <- .Machine$double.eps

            ## Standardized observations that are close to 0 are set to 'eps'.
            if(any(Q < eps, na.rm = TRUE)){
                ## If lambda == 0.5 * dimension, there is another singularity.
                if(abs(lambda - 0.5 * d) < eps){
                    stop("Unhandled singularity: Some standardized observations are close to 0 (< ",
                         eps, ") and lambda is close to 0.5 * dimension!\n")
                }else{
                    Q[Q < eps] <- eps
                    Q[Q == 0] <- eps
                    warning("Singularity: Some standardized observations are close to 0 (< ",
                            sprintf("% .6E", eps),")!\n",
                            "Observations set to ",sprintf("% .6E", eps),".\n",immediate. = TRUE)
                }
            }
            log.top <- .besselM3((lambda - d/2), sqrt((psi + skewness.norm) * (chi + Q)),
                                logvalue = T) + skewness.scaled
            log.bottom <- (d/2 - lambda) * log(sqrt((psi + skewness.norm) * (chi + Q)))
            log.const.top <- d * log(psi)/2 + (1 - lambda) * log(2) +
                (d/2 - lambda) * log(1 + skewness.norm/psi)
            log.const.bottom <- (d/2) * log(2 * pi) + lgamma(lambda) + 0.5 * log(det.sigma)
            out <- log.const.top + log.top - log.const.bottom - log.bottom
        }else{
            out <- NA
        }
    }
    if (!logvalue){
        out <- exp(out)
    }
    return(out)
}
### <---------------------------------------------------------------------->



### <======================================================================>
".is.gaussian" <- function(object)
{
    return(ghyp.name(object, abbr = TRUE) == "Gauss"  ||
           (object@psi == Inf && object@chi == Inf))
}
### <---------------------------------------------------------------------->



### <======================================================================>
".is.student.t" <- function(object, symmetric = NULL)
{
    if(is.null(symmetric)){
        return(ghyp.name(object, abbr = FALSE, skew.attr = FALSE) == "Student-t")
    }else if(symmetric){
        return(ghyp.name(object, abbr = FALSE, skew.attr = TRUE) == "Symmetric Student-t")
    }else if(!symmetric){
        return(ghyp.name(object, abbr = FALSE, skew.attr = TRUE) == "Asymmetric Student-t")
    }else{
        stop("Undefined value received in '.is.student.t'!\n")
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
".is.univariate" <- function(object)
{
    return(object@dimension == 1)
}
### <---------------------------------------------------------------------->



### <======================================================================>
".is.symmetric" <- function(object)
{
    return(all(object@gamma == 0))
}
### <---------------------------------------------------------------------->



### <======================================================================>
".mle.default" <- function(data, pdf, vars, opt.pars = rep(TRUE, length(vars)),
                          transform = NULL, se = FALSE, na.rm = FALSE,
                          silent = FALSE, ...)
{

    if(na.rm){
        data <- data[!is.na(data)]
    }else{
        if(any(is.na(data))){
            stop(".mle.default: The sample contains NAs!\n",
                 "Set na.rm = TRUE to remove the rows which contain NAs.\n")
        }
    }
    ## Sort opt.pars according to vars (-> both vectors must be named)
    opt.pars <- opt.pars[match(names(vars), names(opt.pars))]

    ## Theta contains the parameters intended to be fitted
    theta <- vars[opt.pars]

    ## Theta backup (track parameter difference from optim)
    theta.backup <- theta
    trace.pars <- NULL

    ##<------------   Negative log-Likelihood function adapter ----------->
    negloglik <- function(theta, pdf, tmp.data, transf, const.pars, silent, par.names)
    {
        theta.backup <<- theta

        ## Transformation of the parameters
        for(nam in intersect(names(theta), names(transf))) {
            theta[nam] = do.call(transf[nam], list(theta[nam]))
        }

        pars <- c(theta, const.pars)

        trace.pars <<- rbind(trace.pars, pars[match(par.names, names(pars))])

        pdf.args = c(list(x = tmp.data, logvalue = TRUE), as.list(pars))

        llh <- -sum(do.call(pdf, pdf.args))

        if(!silent){
            print(paste("Llh: ",sprintf("% .14E", -llh), "; Pars: ",
                        paste(sprintf("% .6E", theta), collapse = ", "),
                        sep = ""))
        }
        return(llh)
    }
    ##<------------------------------------------------------------------->

    fit <- try(optim(theta, negloglik, hessian = se, pdf = pdf,
                     tmp.data = data, transf = transform, const.pars = vars[!opt.pars],
                     silent = silent, par.names = names(vars), ...))

    ## 1  indicates that the iteration limit maxit had been reached.
    ## 10 indicates degeneracy of the Nelder-Mead simplex.
    ## 51 indicates a warning from the "L-BFGS-B" method; see '?optim'.
    ## 52 indicates an error from the "L-BFGS-B" method; see '?optim'.

    if(class(fit) == "try-error") {
        convergence <- 100
        hess <- as.numeric(NA)
        n.iter <- as.numeric(NA)
        inv.hess <- matrix(NA)

        message <- fit
        par.ests <- theta.backup
        for(nam in intersect(names(par.ests), names(transform))) {
            par.ests[nam] <- do.call(transform[nam], list(par.ests[nam]))
        }
        vars[opt.pars] <- par.ests

        pdf.args <- c(list(x = data, logvalue = TRUE), as.list(vars))
        ll.max <- -sum(do.call(".dghypuv", pdf.args))
    }else{
        par.ests <- fit$par
        names(par.ests) <- names(theta)
        for(nam in intersect(names(par.ests), names(transform))) {
            par.ests[nam] <- do.call(transform[nam], list(par.ests[nam]))
        }
        vars[opt.pars] <- par.ests
        convergence <- fit$convergence
        n.iter <- fit$counts[1]
        ll.max <- - fit$value
        message <- NULL
        if(se) {
            hess <- fit$hessian
            par.ses <- suppressWarnings(sqrt(diag(hess)))
            inv.hess <- try(solve(hess))
            if(class(inv.hess) == "try-error") {
                warning("Hessian matrix is singular!")
                inv.hess <- matrix(NA, ncol(hess), ncol(hess),
                                   dimnames = dimnames(hess))
            }
            names(par.ses) <- names(par.ests)
            dimnames(hess) <- list(names(par.ests), names(par.ests))
        }else{
            par.ses <- NA
            hess <- NA
            inv.hess <- NA
        }
    }
    list(convergence = convergence, par.ests = vars,
         parameter.variance = inv.hess, ll.max = ll.max, n.iter = n.iter,
         message = message, trace.pars = trace.pars)
}
### <---------------------------------------------------------------------->



### <======================================================================>
".p.default" <- function(q, pdf, pdf.args, lower, upper, ...)
{
    if(missing(upper)){
        int.pars <- list(f = pdf, lower = lower, upper = q)
    }else{
        int.pars <- list(f = pdf, lower = q, upper = upper)
    }
    tmp.prob <- try(do.call("integrate", c(pdf.args, int.pars, list(...))))
    if(class(tmp.prob) == "try-error"){
        warning("Failed to determine probability with 'q = ", q,
                "'!\nMessage: ", as.character(tmp.prob), "\n")
        return(NA)
    }else{
        return(tmp.prob$value)
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
".q.default" <- function(p, pdf, pdf.args, interval, p.lower, ...)
{
    if(p > 0 & p < 1){
        dist.func <- function(x, pdf, p.args, p, p.lower, ...){
            ret.val <- .p.default(q = x, pdf = pdf, pdf.args = p.args,
                                 lower = p.lower, ...)
            return(ret.val - p)
        }
        tmp.quantile <- try(uniroot(dist.func, interval = interval, pdf = pdf,
                                    p.args = pdf.args, p = p,
                                    p.lower = p.lower, ...))

        if(class(tmp.quantile) == "try-error"){
            warning("Failed to determine quantile with 'probs = ", p,
                    "'!\nMessage: ", as.character(tmp.quantile), "\n")
            return(NA)
        }else{
            return(tmp.quantile$root)
        }
    }else if(p == 0){
        return(p.lower)
    }else if(p == 1){
        return(+Inf)
    }else{
        return(NA)
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
".t.transform" <- function(lambda)
{
    return(-1 - exp(lambda))
}
### <---------------------------------------------------------------------->



### <======================================================================>
".inv.t.transform" <- function(lambda.transf)
{
    return(log(-1 - lambda.transf))
}
### <---------------------------------------------------------------------->



### <======================================================================>
.test.ghyp <- function(object, case = c("ghyp", "univariate", "multivariate"))
{
    case = match.arg(case)

    if(!is(object, "ghyp")){
        stop("Object must be of class 'ghyp'!\n")
    }

    if(case == "univariate"){
        if(!.is.univariate(object)){
            stop("'ghyp' object must be univariate (i.e. dimension == 1)!\n")
        }
    }else if(case == "multivariate"){
        if(object@dimension <= 1){
            stop("'ghyp' object must be multivariate (i.e. dimension > 1)!\n")
        }
    }
}
### <---------------------------------------------------------------------->
