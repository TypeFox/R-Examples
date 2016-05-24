### <======================================================================>
"dghyp" <- function(x, object = ghyp(), logvalue = FALSE)
{
    .test.ghyp(object, case = "ghyp")
    x <- .check.data(x, case = if(object@dimension > 1) "mv" else "uv",
                     na.rm = FALSE, fit = FALSE, dim = object@dimension)

    if(.is.univariate(object)){

        d.raw <- rep(NA, length(x))
        d.raw[x == Inf | x == -Inf] <- 0
        d.raw[is.nan(x)] <- NaN

        if(.is.gaussian(object)){
            d.finite <- dnorm(as.vector(x[is.finite(x)]), mean = object@mu,
                              sd = object@sigma, log = logvalue)
        }else{
            d.finite <- unname(.dghypuv(as.vector(x[is.finite(x)]), lambda = object@lambda,
                                        chi = object@chi, psi = object@psi,
                                        mu = object@mu, sigma = object@sigma,
                                        gamma = object@gamma, logvalue = logvalue))
        }
        d.raw[is.finite(x)] <- d.finite

        return(d.raw)
    }else{
        if(.is.gaussian(object)){
            d <- object@dimension
            log.const.bottom <- d / 2 * log(2 * pi) + 0.5 * log(det(vcov(object)))
            log.top <- -0.5 * mahalanobis(x, center = mean(object), cov = vcov(object))
            if(logvalue){
                return(log.top - log.const.bottom)
            }else{
                return(exp(log.top - log.const.bottom))
            }
        }else{
            return(unname(as.vector(.dghypmv(x, lambda = object@lambda,
                                             chi = object@chi, psi = object@psi,
                                             mu = object@mu, sigma = object@sigma,
                                             gamma = object@gamma, logvalue = logvalue))))
        }

    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
"ESghyp" <- function(alpha, object = ghyp(), distr = c("return", "loss"), ...)
{
    distr <- match.arg(distr)

    .test.ghyp(object, case = "univariate")

    alpha <- .check.data(alpha, na.rm = FALSE, fit = FALSE, dim = 1)

    if(.is.gaussian(object)){
        if(distr == "return"){
            return(mean(object) - sqrt(vcov(object)) * dnorm(qnorm(1 - alpha)) / (alpha))
        }else{                          # For losses
            return(mean(object) + sqrt(vcov(object)) * dnorm(qnorm(alpha)) / (1 - alpha))
        }

    }else if(.is.student.t(object, symmetric = TRUE) & object@parametrization == "alpha.bar"){
        nu <- coef(object)$nu
        sigma.t <- sqrt((nu - 2) / nu) * as.vector(object@sigma)
        if(distr == "return"){
            alpha <- 1 - alpha
            return(object@mu - sigma.t * dt(qt(alpha, df = nu), df = nu) /
                   (1 - alpha) * (nu + qt(alpha, df = nu)^2) / (nu - 1))
        }else{
            return(object@mu + sigma.t * dt(qt(alpha, df = nu), df = nu) /
                   (1 - alpha) * (nu + qt(alpha, df = nu)^2) / (nu - 1))
        }
    }else{
        value.raw <- qghyp(alpha, object,...)

        if(all(is.na(value.raw))){
            return(value.raw)
        }

        pdf.args <- list(lambda = object@lambda, chi = object@chi, psi = object@psi,
                         mu = object@mu, sigma = object@sigma, gamma = object@gamma)

        value.es <- matrix(value.raw[!is.na(value.raw)], ncol = 1)

        if(distr == "return"){
            value.es <- apply(value.es, MARGIN = 1, FUN = .p.default,
                              pdf = ".integrate.moment.ghypuv",
                              lower = -Inf, pdf.args = pdf.args)
            value.es <- value.es / alpha
        }else{
            value.es <- apply(value.es, MARGIN = 1, FUN = .p.default,
                              pdf = ".integrate.moment.ghypuv",
                              upper = Inf, pdf.args = pdf.args)
            value.es <- value.es / (1 - alpha)
        }

        value.raw[!is.na(value.raw)] <- value.es

        return(value.es)
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
"ghyp.data" <- function(object)
{
    .test.ghyp(object, case = "ghyp")

    if(length(object@data) != 0){
      if(.is.univariate(object)){
        return(as.vector(object@data))
      }else{
        return(object@data)
      }
    }else{
      return(NULL)
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
"ghyp.dim" <- function(object)
{
    .test.ghyp(object, case = "ghyp")
    return(object@dimension)

}
### <---------------------------------------------------------------------->



### <======================================================================>
"ghyp.fit.info" <- function(object)
{
  if(!is(object, "mle.ghyp")){
    stop("Object is not of class 'mle.ghyp'!")
  }

  if(length(object@trace.pars) == 0){
    trace.pars <- NULL
  }else{
    if(.is.univariate(object)){
      trace.pars <- as.data.frame(object@trace.pars)
    }else{
      trace.pars <- object@trace.pars
    }
  }

  if(.is.univariate(object)){
    return(list(logLikelihood = object@llh,
                aic = object@aic,
                fitted.params = object@fitted.params,
                converged = object@converged,
                n.iter = object@n.iter,
                error.code = object@error.code,
                error.message = object@error.message,
                parameter.variance = object@parameter.variance,
                trace.pars = trace.pars))

  }else{
    return(list(logLikelihood = object@llh,
                aic = object@aic,
                fitted.params = object@fitted.params,
                converged = object@converged,
                n.iter = object@n.iter,
                error.code = object@error.code,
                error.message = object@error.message,
                trace.pars = trace.pars))

  }
}
### <---------------------------------------------------------------------->



### <======================================================================>
"ghyp.moment" <- function(object, order = 3:4,
                          absolute = FALSE, central = TRUE, ...)
{
    .test.ghyp(object, case = "univariate")

    ##     if(.is.gaussian(object)){
    ##         stop("'ghyp.moment' does not work for gaussian objects!")
    ##     }

    if(!is.vector(order)){
        stop("Argument 'order' must be a vector!")
    }

    if(!all(is.finite(order))){
        stop("Argument 'order' must not contain non-finite values!")
    }

    if(any(order %% 1 != 0) & !absolute){
        stop("Argument absolute must be 'TRUE' when order is not integer!")
    }

    if(any(order < 0)){
        stop("'order' must be greater or equal to 0!")
    }

    ## Recursive moment generating function for integer orders
    ## of 'ghyp', 'hyp' and 'NIG' distributions.
    ## This recursion was introduced in a working paper  of D. Scott,
    ## D. Wuertz, and Th. Tran which was presented
    ## in July 2008 at the R/Rmetrics workshop.
    moment.scott <- function(object, order)
    {

        moment.mu.0 <- function(par.list, order)
        {
            a.func <- function(k, l)
            {
                return(factorial(k) /
                       (factorial(k - l) * factorial(2 * l - k) * 2^(k - l)))
            }

            z.func <- function(lambda, x){
                return(.besselM3(lambda, x) / x^lambda)
            }

            eta.0  <- par.list$delta * sqrt(par.list$alpha^2 - par.list$beta^2)
            c.lambda <- 1/z.func(par.list$lambda, eta.0)
            l <- floor((order + 1)/2)
            idx <- l:order
            return(c.lambda * sum(a.func(order, idx) * par.list$delta^(2 * idx) *
                                  par.list$beta^(2 * idx - order) *
                                  z.func(par.list$lambda + idx, eta.0)))
        }

        ab.par.list <- coef(object, type = "alpha.delta")

        tmp.moment <- 0
        for(i in 0:order){
            tmp.addend <- choose(order, i) * ab.par.list$mu^(order - i) *
                moment.mu.0(ab.par.list, i)
            tmp.moment <- tmp.moment + tmp.addend
        }
        return(tmp.moment)
    }

    internal.moment <- function(x, lambda, chi, psi, mu, sigma, gamma,
                                tmp.order, location.offset, transf.fun){
        return(.dghypuv(x,lambda = lambda, chi = chi, psi = psi, mu = mu,
                        sigma = sigma, gamma = gamma) * transf.fun(x - location.offset)^tmp.order)

    }

    internal.moment.gaussian <- function(x, tmp.order, mu, sigma, location.offset, transf.fun){
        return(transf.fun(x - location.offset)^tmp.order * dnorm(x, mu, sigma))
    }

    if(central){
        location.offset <- mean(object)
    }else{
        location.offset <- 0
    }

    if(absolute){
        transf.fun <- abs
    }else{
        transf.fun <- identity
    }

    result <- numeric(length(order))

    for (i in 1:length(order))
    {
        if(.is.gaussian(object)){
            if(absolute | !central){
                tmp.result <- try(integrate(internal.moment.gaussian, lower = -Inf, upper = Inf,
                                            tmp.order = order[i], mu = mean(object), sigma = sqrt(vcov(object)),
                                            transf.fun = transf.fun, location.offset = location.offset, ...)$value)
                if(class(tmp.result) == "try-error"){
                    result[i] <- NA
                }else{
                    result[i] <- tmp.result
                }
            }else{    # non-absolute central moments of integer order:
                if(order[i] == 0){
                    ## 'order[i] %% 1 == 0' implies that non-absolute
                    ## moments are requested
                    result[i] <- 1
                }else if(central){
                    if(order[i] %% 2 == 0){ # even moments
                        k <- order[i] / 2
                        result[i] <- factorial(order[i]) / (2^k * factorial(k))
                    }else{              # odd moments
                        result[i] <- 0
                    }
                }
            }
        }else{
            if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) %in%
               c("ghyp", "hyp", "NIG") & (order[i] %% 1 == 0) & !absolute & !central)
            {                           # analytical formulas exist
                result[i] <- moment.scott(object, order[i])
            }else{
                                        # do it numerically
                tmp.result <- try(integrate(internal.moment, -Inf, Inf, tmp.order = order[i],
                                            lambda = object@lambda, chi = object@chi,
                                            psi = object@psi, mu = object@mu,
                                            sigma = object@sigma, gamma = object@gamma,
                                            transf.fun = transf.fun,
                                            location.offset = location.offset, ...)$value)
                if(class(tmp.result) == "try-error"){
                    result[i] <- NA
                }else{
                    result[i] <- tmp.result
                }
            }
        }
    }
    return(result)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"ghyp.name" <- function(object, abbr = FALSE, skew.attr = TRUE)
{
    .test.ghyp(object, case = "ghyp")
    if(!abbr && skew.attr){
      return(object@model[1])
    }else if(abbr && skew.attr){
      return(object@model[3])
    }else if(abbr && !skew.attr){
      return(object@model[4])
    }else{
      return(object@model[2])
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
"ghyp.omega" <- function(L, object = ghyp(), ...)
{
    ## Only univariate 'ghyp' objects are allowed
    .test.ghyp(object, case = "univariate")

    call.part <- function(x, L, obj)
        return((x - L) * dghyp(x, obj))

    put.part <- function(x, L, obj)
        return((L - x) * dghyp(x, obj))

    int.num <- numeric(length(L))
    int.denom <- numeric(length(L))
    for(i in 1:length(L)){
        ##         int.num[i] <- integrate(pghyp, object = object, lower = L[i],
        ##                                 upper = Inf, lower.tail = FALSE, ...)$value
        ##         int.denom[i] <- integrate(pghyp, object = object, lower = -Inf,
        ##                                   upper = L[i], ...)$value
        tmp.num <- try(integrate(call.part, lower = L[i], upper = Inf,
                                 obj = object, L = L[i], ...))
        if(class(tmp.num) == "try-error"){
            warning("Integral 'int_L^\\infty (x - L) * dghyp(x, obj) dx' did ",
                    "not converge for L =", L[i], "!")
            int.num[i] <- NA
        }else{
            int.num[i] <- tmp.num$value
        }

        tmp.denom <- try(integrate(put.part, lower = -Inf, upper = L[i],
                                   obj = object, L = L[i], ...))
        if(class(tmp.denom) == "try-error"){
            warning("Integral 'int_-\\infty^L (L - x) * dghyp(x, obj) dx' ",
                    "did not converge for L =", L[i], "!")
            int.denom[i] <- NA
        }else{
            int.denom[i] <- tmp.denom$value
        }
    }
    return(int.num / int.denom)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"ghyp.kurtosis" <- function(object)
{
    ## Only univariate 'ghyp' objects are allowed
    .test.ghyp(object, case = "univariate")

    if(.is.gaussian(object))
        return(3)

    second.third.fourth.moment <- ghyp.moment(object, order = 2:4, absolute = FALSE,
                                              central = FALSE)

    central.fourth.moment <- second.third.fourth.moment[3] -
        4 * mean(object) * second.third.fourth.moment[2] +
            6 * mean(object)^2 * second.third.fourth.moment[1] -
                3 * mean(object)^4

    return(central.fourth.moment / vcov(object)^2)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"lik.ratio.test" <- function(x, x.subclass, conf.level = 0.95)
{
    if(!is(x, "mle.ghyp") | !is(x.subclass, "mle.ghyp")){
        stop("Objects are not of class 'mle.ghyp'!")
    }

    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                                 conf.level < 0 || conf.level > 1)){
        stop("'conf.level' must be a single number between 0 and 1")
    }

    df <- sum(x@fitted.params) - sum(x.subclass@fitted.params)

    if(df < 1){
        stop("Degree of freedom equals zero. Consult the help for more information!")
    }

    L.stat.chisq <- 2 * (logLik(x) - logLik(x.subclass))

    chisq.quantile <- qchisq(conf.level, df = df)

    p.value <- pchisq(L.stat.chisq, df = df, lower.tail = FALSE) # P(X > L.stat.chisq)

    L.stat <- exp(logLik(x.subclass)-logLik(x))

    return(list(statistic = c(L = L.stat), p.value = p.value, df = df,
                H0 = L.stat.chisq <= chisq.quantile))
}
### <---------------------------------------------------------------------->



### <======================================================================>
"pghyp" <- function(q, object = ghyp(), n.sim = 10000, subdivisions = 200,
                    rel.tol = .Machine$double.eps^0.5, abs.tol = rel.tol,
                    lower.tail = TRUE)
{
    .test.ghyp(object, case = "ghyp")
    q.raw  <- .check.data(q, case = if(object@dimension > 1) "mv" else "uv",
                          na.rm = FALSE, fit = FALSE, dim = 1)

    if(.is.univariate(object)){
        if(.is.gaussian(object)){
            return(pnorm(q, mean = object@mu, sd = as.vector(object@sigma),
                         lower.tail = lower.tail))
        }else if(.is.student.t(object, symmetric = TRUE) & object@parametrization == "alpha.bar"){
            nu <- coef(object)$nu
            return(pt((q - object@mu) / (sqrt((nu - 2) / nu) * as.vector(object@sigma)),
                      df = nu, lower.tail = lower.tail))
        }else{
            q.finite <- q.raw[is.finite(q.raw)]
            q.mat <- matrix(q.finite, ncol = 1)

            p.raw <- rep(NA, length(q.raw))

            pdf.args <- list(lambda = object@lambda, chi = object@chi, psi = object@psi,
                             mu = object@mu, sigma = object@sigma, gamma = object@gamma)


            if(lower.tail){
                p.raw[q.raw == -Inf] <- 0
                p.raw[q.raw == Inf] <- 1
                value <- apply(q.mat, MARGIN = 1, FUN = .p.default, pdf = ".dghypuv",
                               lower = -Inf, pdf.args = pdf.args, subdivisions = subdivisions,
                               rel.tol = rel.tol, abs.tol = abs.tol)
            }else{
                p.raw[q.raw == -Inf] <- 1
                p.raw[q.raw == Inf] <- 0
                value <- apply(q.mat, MARGIN = 1, FUN = .p.default, pdf = ".dghypuv",
                               upper = Inf, pdf.args = pdf.args, subdivisions = subdivisions,
                               rel.tol = rel.tol, abs.tol = abs.tol)
            }

            p.raw[is.finite(q.raw)] <- value
            return(as.vector(p.raw))
        }
    }else{
        sim.data <- rghyp(n.sim, object)
        if(lower.tail){
            compare.fun <- '<'
        }else{
            compare.fun <- '>'
        }

        eval.comparison <- function(q.raw, sim.data, n.sim, tmp.compare.fun)
        {
            return(sum(apply(apply(sim.data, MARGIN = 1, FUN = tmp.compare.fun, q.raw),
                             MARGIN = 2, FUN = all)) / n.sim)

        }

        return(apply(q.raw, MARGIN = 1, FUN = eval.comparison,
                     tmp.compare.fun = compare.fun, sim.data = sim.data, n.sim = n.sim))
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
"portfolio.optimize" <- function(object,
                                 risk.measure = c("sd", "value.at.risk", "expected.shortfall"),
                                 type = c("minimum.risk", "tangency", "target.return"),
                                 level = 0.95, distr = c("loss", "return"),
                                 target.return = NULL, risk.free = NULL,
                                 silent = FALSE, ...)
{
    .test.ghyp(object, case = "multivariate")

    risk.measure <- match.arg(risk.measure)
    type <- match.arg(type)
    distr <- match.arg(distr)

    if(type == "target.return" & is.null(target.return)){
        stop("If optimization problem is of type 'target.return' the ",
             "argeument 'target.return' must be specified!")

    }
    if(type == "tangency" & is.null(risk.free)){
        stop("If optimization problem is of type 'tangency' the ",
             "argeument 'risk.free' must be specified!")

    }

    if(type == "tangency" & !.is.symmetric(object) & risk.measure != "sd"){
        stop("Type 'tangency' optimization problem ",
             "not implemented for non-symmetric distributions and ",
             "risk measure not equal to 'sd'!")

    }

    ## This function handles all optimization problems if either
    ## risk is measured in terms of variance or if the distribution
    ## is symmetric.
    mean.variance <- function(means, variance, r.f = NULL, target.return = NA,
                              type = c("tangency", "minimum.risk", "target.return"))
    {
        type <- match.arg(type)
        if(type == "tangency"){
            weights.unnormed <- as.vector(solve(variance) %*% (means - r.f))
            weights <- weights.unnormed / sum(weights.unnormed)
        }else{
            dimension <- nrow(variance)
            lin.sys <- cbind(variance, rep(-1, dimension))
            lin.sys <- rbind(lin.sys, c(rep(1, dimension), 0))

            if(type == "minimum.risk"){
                weights.lambda <- solve(lin.sys, c(rep(0, dimension), 1))
                weights <- weights.lambda[1:dimension]
            }else{                      # target.return
                lin.sys <- cbind(lin.sys, c(- means, 0))
                lin.sys <- rbind(lin.sys, c(means, 0, 0))
                weights.lambda.1.2 <- solve(lin.sys, c(rep(0, dimension),
                                                       1, target.return))
                weights <- weights.lambda.1.2[1:dimension]
            }

        }
        ptf.mean <- sum(weights * means)
        ptf.sd <- as.numeric(sqrt(weights %*% variance %*% weights))
        return(list(ptf.mean = ptf.mean, ptf.sd = ptf.sd, weights = weights))
    }
    ## End of function 'mean.variance'.


    if(distr == "return"){
        object <- transform(object, multiplier = diag(-1, ghyp.dim(object)))
        if(!is.null(target.return)){
            target.return <- - 1 * target.return
        }
    }

    if(risk.measure == "sd" | .is.symmetric(object))
    {
        mv.opt <- mean.variance(mean(object), vcov(object),
                                type = type, r.f = risk.free,
                                target.return = target.return)

        opt.weights <- mv.opt$weights

        ptf.dist <- transform(object, multiplier = t(opt.weights))

        risk <- switch(risk.measure,
                       sd = sqrt(vcov(ptf.dist)),
                       value.at.risk = qghyp(level, ptf.dist),
                       expected.shortfall = ESghyp(level, ptf.dist, distr = "loss"))
        converged <- TRUE

        n.iter <- 0
        message <- ""
    }else{
        ## Risk measure != "sd" and the distribution object is non-symmetric.
        if(type == "target.return")
        {
            objective.target.return <- function(pars, object, target.return, level,
                                                risk.measure, silent = TRUE)
            {
                weight.2 <- (target.return + mean(object)[1] * (sum(pars) - 1) -
                             sum(pars * mean(object)[3:ghyp.dim(object)])) /
                                 (mean(object)[2] - mean(object)[1])

                weight.1 <- 1 - sum(pars) - weight.2
                weights <- c(weight.1, weight.2, pars)

                backup.weights <<- weights

                ptf.dist <- transform(object, multiplier = t(weights))

                if(risk.measure == "value.at.risk"){
                    risk <- qghyp(level, ptf.dist)
                }else if(risk.measure == "expected.shortfall"){
                    risk <- ESghyp(level, ptf.dist, distr = "loss")
                }else{
                    stop("Unknown risk measure: ", risk.measure)
                }

                if(!silent){
                    print(paste("risk: ",sprintf("% .7E", risk),
                                "; weights: ",
                                paste(sprintf("% .4E", weights), collapse = ", "),
                                sep = ""))
                }

                return(risk)
            }

            initial.weights <- rep(1 / ghyp.dim(object), ghyp.dim(object) - 2)
            backup.weights <- initial.weights

            opt.ptf <- try(optim(initial.weights, objective.target.return,
                                 object = object,
                                 target.return = target.return,
                                 level = level, risk.measure = risk.measure,
                                 silent = silent, ...))

            if(class(opt.ptf) == "try-error")
            {
                converged <- FALSE
                risk <- NA
                opt.weights <- backup.weights
                ptf.dist <- transform(object, multiplier = t(opt.weights))
                message <- opt.ptf
            }else{
                weight.2 <- (target.return + mean(object)[1] * (sum(opt.ptf$par) - 1) -
                             sum(opt.ptf$par * mean(object)[3:ghyp.dim(object)])) /
                                 (mean(object)[2] - mean(object)[1])

                weight.1 <- 1 - sum(opt.ptf$par) - weight.2

                opt.weights <- c(weight.1, weight.2, opt.ptf$par)

                ptf.dist <- transform(object, multiplier = t(opt.weights))
                risk <- switch(risk.measure,
                               value.at.risk = qghyp(level, ptf.dist),
                               expected.shortfall = ESghyp(level, ptf.dist, distr = "loss"))

                if(opt.ptf$convergence == 0){
                    converged <- TRUE
                }else{
                    converged <- FALSE
                }
                message <- opt.ptf$message
            }


            n.iter <- opt.ptf$counts[1]

        }else{
###<------------- type == "minimum.risk" or "tangency" --------------->
            objective <- function(pars, object, type, r.f, level,
                                  risk.measure, silent = TRUE)
            {
                weights <- c(1 - sum(pars), pars)

                backup.weights <<- weights
                ptf.dist <- transform(object, multiplier = t(weights))

                if(risk.measure == "value.at.risk"){
                    risk <- qghyp(level, ptf.dist)
                }else if(risk.measure == "expected.shortfall"){
                    risk <- ESghyp(level, ptf.dist, distr = "loss")
                }else{
                    stop("Unknown risk measure: ", risk.measure)
                }

                if(type == "tangency"){
                    func.val <- - (mean(ptf.dist) - r.f) / risk
                }else if(type == "minimum.risk"){
                    func.val <- risk
                }else{
                    stop("Unknown optimization problem: ", type)
                }
                if(!silent){
                    print(paste("fct value: ",sprintf("% .7E", func.val),
                                "; weights: ",
                                paste(sprintf("% .4E", weights), collapse = ", "),
                                sep = ""))
                }
                return(func.val)
            }

            initial.weights <- rep(1 / ghyp.dim(object), ghyp.dim(object) - 1)
            backup.weights <- initial.weights

            opt.ptf <- try(optim(initial.weights, objective, object = object, type = type,
                                 r.f = risk.free, level = level, risk.measure = risk.measure,
                                 silent = silent, ...))
            if(class(opt.ptf) == "try-error")
            {
                converged <- FALSE
                risk <- NA
                opt.weights <- backup.weights
                ptf.dist <- transform(object, multiplier = t(opt.weights))
                message <- opt.ptf
            }else{

                opt.weights <- c(1 - sum(opt.ptf$par), opt.ptf$par)

                ptf.dist <- transform(object, multiplier = t(opt.weights))
                risk <- switch(risk.measure,
                               value.at.risk = qghyp(level, ptf.dist),
                               expected.shortfall = ESghyp(level, ptf.dist, distr = "loss"))

                if(opt.ptf$convergence == 0){
                    converged <- TRUE
                }else{
                    converged <- FALSE
                }
                message <- opt.ptf$message
            }

            n.iter <- opt.ptf$counts[1]

        }
    }

    if(distr == "return"){
        ptf.dist <- transform(ptf.dist, multiplier = -1)
    }

    return(list(portfolio.dist = ptf.dist, risk.measure = risk.measure,
                risk = risk, opt.weights = unname(opt.weights),
                converged = converged, message = message,
                n.iter = n.iter))
}
### <---------------------------------------------------------------------->



### <======================================================================>
"qghyp" <- function(p, object = ghyp(), method = c("integration", "splines"),
                    spline.points = 200, subdivisions = 200,
                    root.tol = .Machine$double.eps^0.5,
                    rel.tol = root.tol^1.5, abs.tol = rel.tol)

{
    ## Only univariate 'ghyp' objects are allowed
    .test.ghyp(object, case = "univariate")

    if(.is.gaussian(object)){
        return(qnorm(p, mean = object@mu, sd = object@sigma))
    }else if(.is.student.t(object, symmetric = TRUE) & object@parametrization == "alpha.bar"){
        nu <- coef(object)$nu
        return(qt(p, df = nu) * sqrt((nu - 2) / nu) * as.numeric(object@sigma) + object@mu)
    }

    p.raw  <- .check.data(p, na.rm = FALSE, fit = FALSE, dim = 1)

    method <- match.arg(method)

    p.raw[p.raw < 0 | p.raw > 1] <- NaN
    p.raw[p.raw == 1] <- Inf
    p.raw[p.raw == 0] <- -Inf

    ## If only non-finite quantiles are passed return NA, NaN, Inf, -Inf
    p <- p.raw[is.finite(p.raw)]
    if(length(p) == 0){
        return(p.raw)
    }

    ##<----   Use Newton's method to find the range of the quantiles ---->
    internal.bisection <- function(object, p, tol, rel.tol, abs.tol, subdivisions)
    {
        iter <- 0
        range.found <- FALSE

        step.size <- sqrt(vcov(object))
        if(!is.finite(step.size)){
            step.size <- coef(object, type = "chi.psi")$sigma / 2
        }

        q.0 <- mean(object)
        if(!is.finite(q.0)){
            q.0 <- coef(object, type = "chi.psi")$mu
        }

        q.range <- c(q.0 - step.size, q.0 + step.size)

        while(!range.found & iter < 100){
            iter <- iter + 1

            p.range <- pghyp(q = q.range, object, rel.tol = rel.tol, abs.tol = abs.tol,
                             subdivisions = subdivisions) - p

            if(any(is.na(p.range))){
                warning("Unable to determine interval where the quantiles are in-between.\n",
                        "Perhaps the skewness is too large!")
                return(NA)
            }

            lower <- p.range[1]
            upper <- p.range[2]

            ##      cat("lower: ", lower,";  upper : ", upper, "\n")
            if(upper < 0 & lower < 0){
                q.range[1] <- q.range[2]
                q.range[2] <- q.range[2] + step.size
                next
            }
            if(upper > 0 & lower > 0){
                q.range[2] <- q.range[1]
                q.range[1] <- q.range[1] - step.size
                next
            }
            if(upper > 0 & lower < 0){
                range.found <- TRUE
            }
        }
        if(iter >= 100){
            warning("Unable to determine interval where the quantiles are in-between.\n",
                    "Perhaps the skewness is too large!")
        }

        q.root <- .q.default(p, pdf = ".dghypuv", pdf.args = coef(object, type = "chi.psi"),
                             interval = q.range, tol = root.tol,
                             p.lower = -Inf, rel.tol = rel.tol, abs.tol = abs.tol,
                             subdivisions = subdivisions)
        return(q.root)
    }
    ##<---------------- end of Newton iteration  ---------------------->

    if(length(p) == 1){
        ## If a single quantile is requested use the newton method anyway
        value <- internal.bisection(object, p, root.tol, rel.tol, abs.tol, subdivisions)
        p.raw[is.finite(p.raw)] <- as.numeric(value)
        return(p.raw)
    }else if(length(p) == 2){
        ## If two quantiles are requested use the newton method anyway
        value1 <- internal.bisection(object, p[1], root.tol, rel.tol, abs.tol, subdivisions)
        value2 <- internal.bisection(object, p[2], root.tol, rel.tol, abs.tol, subdivisions)
        p.raw[is.finite(p.raw)] <- c(value1, value2)
        return(p.raw)
    }else{
        ## If more than two quantiles are requested use the newton method
        ## to find the range where the quantiles can be found.
        q.min <- internal.bisection(object, min(p), root.tol, rel.tol, abs.tol, subdivisions)
        q.max <- internal.bisection(object, max(p), root.tol, rel.tol, abs.tol, subdivisions)

        interval <- c(q.min, q.max)

        if(any(is.na(interval))){ # -> Failed to determine bounds for the quantiles
            p.raw[is.finite(p.raw)] <- NA
            return(p.raw)
        }

        ## Extend the interval by 10 percent so that 'uniroot' does not crash
        interval <- c(interval[1] - 0.01 * diff(range(interval)),
                      interval[2] + 0.01 * diff(range(interval)))

        if(method == "integration"){    # Integration method

            pdf.args <- coef(object, type = "chi.psi")
            p <- matrix(p, ncol = 1)
            value <- apply(p, MARGIN = 1, FUN = .q.default, pdf = ".dghypuv",
                           pdf.args = pdf.args, interval = interval, tol = root.tol,
                           p.lower = -Inf, rel.tol = rel.tol, abs.tol = abs.tol,
                           subdivisions = subdivisions)
        }else{                          # Splines method
            interval.seq <- seq(min(interval), max(interval), length = spline.points)
            ## Compute the distribution function to be interpolated by splines
            p.interval <- pghyp(q = interval.seq, object, rel.tol = rel.tol,
                                abs.tol = abs.tol, subdivisions = subdivisions)

            ## Spline function
            spline.distribution.func <- splinefun(interval.seq, p.interval)

            ## root function:   condition: quantile.root.func == 0
            quantile.root.func <- function(x, tmp.p){
                spline.distribution.func(x) - tmp.p
            }

            value <- p

            for(i in 1:length(p)){
                value[i] <- uniroot(quantile.root.func, interval = interval,
                                    tol = root.tol, tmp.p = p[i])$root
            }
        }
        p.raw[is.finite(p.raw)] <- value
        return(p.raw)
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
"qqghyp" <- function(object, data = ghyp.data(object),
                     gaussian = TRUE, line = TRUE,
                     main = "Generalized Hyperbolic Q-Q Plot",
                     xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                     ghyp.pch = 1, gauss.pch = 6, ghyp.lty = "solid",
                     gauss.lty = "dashed", ghyp.col = "black",
                     gauss.col = "black",
                     plot.legend = TRUE, location = "topleft", legend.cex = 0.8,
                     spline.points = 150, root.tol = .Machine$double.eps^0.5,
                     rel.tol = root.tol, abs.tol = root.tol^1.5, add = FALSE, ...)
{
    .test.ghyp(object, case = "univariate")

    data <- .check.data(data, case = "uv", na.rm = T, fit = TRUE, dim = 1)

    ## compute quantiles
    ghyp.q <- qghyp(ppoints(length(data)), object, method = "splines",
                    spline.points = spline.points, root.tol = root.tol,
                    rel.tol = rel.tol, abs.tol = abs.tol)[order(order(data))]

    ## plot ghyp quantiles
    if(add){
        points(ghyp.q, data, pch = ghyp.pch, col = ghyp.col, ...)
    }else{
        plot(ghyp.q, data, xlab = xlab, ylab = ylab, pch = ghyp.pch,
             col = ghyp.col, main = main, ...)
    }

    if(gaussian){
        gauss.q <- qnorm(ppoints(length(data)), mean = mean(data),
                         sd = sd(data))[order(order(data))]
        points(gauss.q, data, pch = gauss.pch, col = gauss.col)
    }

    if(line){
        abline(lm(data ~ ghyp.q), lty = ghyp.lty, col = ghyp.col)
        if(gaussian){
            abline(lm(data ~ gauss.q), lty = gauss.lty, col = gauss.col)
        }
    }

    if(plot.legend && !add){
        if(gaussian){
            legend(location,
                   legend = c(ghyp.name(object, abbr = TRUE, skew.attr = TRUE), "Gaussian"),
                   col = c(ghyp.col, gauss.col), lty = c(ghyp.lty, gauss.lty),
                   cex = legend.cex, pch = c(ghyp.pch, gauss.pch))
        }else{
            legend(location, legend = ghyp.name(object, abbr = TRUE, skew.attr = TRUE),
                   col = ghyp.col, lty = ghyp.lty, cex = legend.cex, pch = ghyp.pch)
        }

    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
"rghyp" <- function(n, object = ghyp())
{
    .test.ghyp(object, case = "ghyp")
    if(.is.univariate(object)){
        if(.is.gaussian(object)){
            return(rnorm(n, mean = object@mu, sd = object@sigma))
        }else if(.is.student.t(object, symmetric = TRUE) & object@parametrization == "alpha.bar"){
            nu <- coef(object)$nu
            return(rt(n, df = nu) * sqrt((nu - 2) / nu) * object@sigma + object@mu)
        }else{
            W <- rgig(n, object@lambda, object@chi, object@psi)
            return(return(object@mu + W * object@gamma + sqrt(W) * object@sigma * rnorm(n)))
        }
    }else{
        Z <- matrix(rnorm(n * object@dimension), ncol = object@dimension)
        A <- chol(object@sigma, pivot = FALSE)
        if(.is.gaussian(object)){
            return(Z %*% A  +  matrix(rep(object@mu, n), ncol = object@dimension, byrow = TRUE))
        }else{
            W <- rgig(n, object@lambda, object@chi, object@psi)
            return(sqrt(W) * (Z %*% A)  +
                   matrix(rep(object@mu, n), ncol = object@dimension, byrow = TRUE) +
                   outer(W, object@gamma))
        }
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
"ghyp.skewness" <- function(object)
{
    ## Only univariate 'ghyp' objects are allowed
    .test.ghyp(object, case = "univariate")

    ## For symmetric objects directly return 0
    if(object@gamma == 0)
        return(0)

    second.third.moment <- ghyp.moment(object, order = 2:3, absolute = FALSE,
                                       central = FALSE)

    central.third.moment <- second.third.moment[2] - 3 * mean(object) *
        second.third.moment[1] + 2 * mean(object)^3

    return(central.third.moment / vcov(object)^(1.5))
}
### <---------------------------------------------------------------------->
