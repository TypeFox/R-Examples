### <======================================================================>
"AIC.mle.ghyp" <- function(object, ..., k = 2)
{
    ghyp.models <- list(...)
    ## Internal function; test whether all objects are of type mle.ghyp
    test.class.mle.ghyp <- function(x){
        if(!is(x, "mle.ghyp")){
            stop("Some of the objects are not of class 'mle.ghyp'!")
        }
    }

    ## Internal function; extract the number of fitted parameters
    nbr.fitted.params <- function(tmp.object){
        tmp <- ghyp.fit.info(tmp.object)
        opt.pars <- tmp$fitted.params
        if(.is.univariate(tmp.object)){
            return(unname(sum(opt.pars)))
        }else{
            if(.is.gaussian(tmp.object)){
                return(unname(tmp.object@dimension +
                              tmp.object@dimension/2 * (tmp.object@dimension + 1) *
                              opt.pars[c("sigma")]))
            }else{
                return(unname(sum(opt.pars[c("alpha.bar","lambda")]) +
                              tmp.object@dimension * sum(opt.pars[c("mu","gamma")]) +
                              tmp.object@dimension/2 * (tmp.object@dimension + 1) *
                              opt.pars[c("sigma")]))
            }
        }
    }

    sapply(ghyp.models, test.class.mle.ghyp)

    ghyp.models <- c(object, ghyp.models)

    nbr.fitted <- sapply(ghyp.models, nbr.fitted.params)

    tmp.aic <- -2 * logLik(object,...) + k * nbr.fitted

    return(tmp.aic)
}
### <---------------------------------------------------------------------->
setMethod("AIC", signature(object = "mle.ghyp"), AIC.mle.ghyp)
### <---------------------------------------------------------------------->


### <======================================================================>
"coef.ghyp" <- function(object, type = c("chi.psi", "alpha.bar", "alpha.delta"))
{

    if(.is.univariate(object)){
        sigma <- as.vector(object@sigma)
    }else{
        sigma <- object@sigma
    }

    if(.is.gaussian(object)){
        return(list(mu = object@mu, sigma = sigma))
    }

    if(missing(type)){
        type <- object@parametrization
    }else{
        type <- match.arg(type)
    }

    if(type == "chi.psi"){
        param.list <- list(lambda = object@lambda,
                           chi = object@chi,
                           psi = object@psi,
                           mu = object@mu,
                           sigma = sigma,
                           gamma = object@gamma)
    }else if(type == "alpha.bar"){
        if(.is.student.t(object) & object@parametrization == "alpha.bar"){
            param.list <- list(lambda = object@lambda,
                               nu = -2 * object@lambda,
                               alpha.bar = object@alpha.bar,
                               mu = object@mu,
                               sigma = sigma,
                               gamma = object@gamma)
        }else{
            if(.is.student.t(object) & object@lambda >= -1){
                stop("Transformation to 'alpha.bar' parametrization not possible in case of Student-t distribution with",
                     " a degree of freedom <= 2!")
            }
            k <- Egig(object@lambda, object@chi, object@psi, func = "x")
            if(.is.univariate(object)){
                transf.sigma <- sqrt(k) * sigma
            }else{
                transf.sigma <- k * sigma
            }
            param.list <- list(lambda = object@lambda,
                               alpha.bar = sqrt(object@chi * object@psi),
                               mu = object@mu,
                               sigma = transf.sigma,
                               gamma = k * object@gamma)
        }
    }else if(type == "alpha.delta"){
        if(.is.univariate(object)){     #Univariate
            sigma <- as.numeric(object@sigma)
            alpha <- sqrt(1/sigma^2 * (object@psi + object@gamma^2/sigma^2))
            beta <- object@gamma / sigma^2
            delta <- sqrt(object@chi * sigma^2)
            param.list <- list(lambda = object@lambda, alpha = alpha, delta = delta, beta = beta,
                               mu = object@mu)
        }else{                          # Multivariate
            dimension <- object@dimension
            inv.sigma <- solve(object@sigma)
            det.sigma <- det(object@sigma)
            alpha <- sqrt(as.numeric(det.sigma^(-1/dimension) *
                                     (object@psi + object@gamma %*% inv.sigma %*% object@gamma)))
            beta <- as.numeric(inv.sigma %*% object@gamma)
            delta <- sqrt(object@chi * det.sigma^(1/dimension))
            Delta <- det.sigma^(-1/dimension) * object@sigma
            param.list <- list(lambda = object@lambda, alpha = alpha, delta = delta, beta = beta,
                               mu = object@mu, Delta = Delta)
        }
    }
    return(param.list)
}
### <---------------------------------------------------------------------->
setMethod("coef", signature(object = "ghyp"), coef.ghyp)
### <---------------------------------------------------------------------->
setMethod("coefficients", signature(object = "ghyp"), coef.ghyp)
### <---------------------------------------------------------------------->


### <======================================================================>
"[.ghyp" <- function(x, i = c(1, 2))
{
    .test.ghyp(x, case = "multivariate")
    i <- as.integer(i)
    if(min(i) < 1 | max(i) > ghyp.dim(x)){
        stop("Dimension mismatch. ghyp dimension = ", ghyp.dim(x),
             "; Input = ",paste(i, collapse=", "), "!\n", sep = "")
    }

    if(length(ghyp.data(x)) == 0){
        data <- NULL
    }else{
        data <- ghyp.data(x)[, i]
    }

    if(length(i)==1){
       sigma <- sqrt(x@sigma[i, i])
    }else{
       sigma <- x@sigma[i, i]
    }

    if(x@parametrization == "alpha.bar"){
        if(.is.student.t(x)){
            return(student.t(nu = -2 * x@lambda, mu = x@mu[i],
                             sigma = sigma, gamma = x@gamma[i],
                             data = data))
        }else{
            return(ghyp(lambda = x@lambda, alpha.bar = x@alpha.bar,
                        mu = x@mu[i], sigma = sigma, gamma = x@gamma[i],
                        data = data))
        }
    }else if(x@parametrization == "chi.psi"){
        return(ghyp(lambda = x@lambda, chi = x@chi, psi = x@psi,
                    mu = x@mu[i], sigma = sigma, gamma = x@gamma[i],
                    data = data))
    }else if(x@parametrization == "alpha.delta"){
        tmp.obj <- ghyp(lambda = x@lambda, chi = x@chi, psi = x@psi,
                        mu = x@mu[i], sigma = sigma, gamma = x@gamma[i])
        return(do.call("ghyp.ad", c(coef(tmp.obj, type = "alpha.delta"), list(data = data))))
    }else if(x@parametrization == "Gaussian"){
        return(gauss(mu = x@mu[i], sigma = sigma, data = data))
    }else{
        stop("Unknown parametrization: ", x@parametrization)
    }


}
### <---------------------------------------------------------------------->
setMethod("[", signature(x = "ghyp", i = "numeric", j = "missing", drop = "missing"), `[.ghyp`)
### <---------------------------------------------------------------------->


### <======================================================================>
"hist.ghyp"   <- function(x, data = ghyp.data(x),
                            gaussian = TRUE, log.hist = F, ylim = NULL,
                            ghyp.col = 1, ghyp.lwd = 1, ghyp.lty = "solid",
                            col = 1, nclass = 30, plot.legend = TRUE,
                            location = if(log.hist) "bottom" else "topright",
                            legend.cex = 1, ...)
{
    if(missing(data) & is.null(ghyp.data(x))){
        stop("'data' must be provided when object 'x' contains no data")
    }

    .test.ghyp(x, "univariate")

    data <- .check.data(data, case = "uv", na.rm = T, fit = TRUE, dim = 1)

    x.gh <- seq(min(data), max(data), length = 2000)
    tmp.d.ghyp <- dghyp(x.gh, x)

    if(is.null(ylim)){
        ylim <- c(0, max(tmp.d.ghyp))
    }

    if(log.hist == TRUE){
        tmp.hist <- suppressWarnings(hist(data, probability = T, plot = F, nclass = nclass, ...))
        ghyp.data <- tmp.hist$breaks[-1] - diff(tmp.hist$breaks)[1]/2
        Density <- tmp.hist$density
        plot(ghyp.data, log(Density), col = col, ...)
        lines(x.gh, log(tmp.d.ghyp), col = ghyp.col, lwd = ghyp.lwd, lty = ghyp.lty)
        if(gaussian == TRUE){
            lines(x.gh, log(dnorm(x.gh, mean = mean(data), sd = sd(data))), col = col, lty = "dashed")
        }
    }else{
        hist(data, ylim = ylim, probability = T, nclass = nclass, ...)
        lines(x.gh, tmp.d.ghyp, col = ghyp.col, lwd = ghyp.lwd, lty = ghyp.lty)
        if(gaussian == TRUE){
            lines(x.gh, dnorm(x.gh, mean = mean(data), sd = sd(data)), col = col, lty = "dashed")
        }
    }
    if(plot.legend == TRUE){

        if(log.hist == TRUE){
            if(gaussian == TRUE){
                tmp.text <- c("Histogramm", ghyp.name(x, abbr = TRUE, skew.attr = TRUE), "Gaussian")
                tmp.col <- c(col, ghyp.col, col)
                tmp.lty <- c(NA, ghyp.lty, "dashed")
                tmp.pch <- c(1, NA, NA)
            }else{
                tmp.text <- c("Histogramm", ghyp.name(x, abbr = TRUE, skew.attr = TRUE))
                tmp.col <- c(col, ghyp.col)
                tmp.lty <- c(NA, ghyp.lty)
                tmp.pch <- c(1, NA)
            }
            legend(location, legend = tmp.text, col = tmp.col,
                   lty = tmp.lty, pch = tmp.pch, cex = legend.cex)
        }else{
            if(gaussian == TRUE){
                tmp.text <- c(ghyp.name(x, abbr = TRUE, skew.attr = TRUE), "Gaussian")
                tmp.col <- c(ghyp.col, col)
                tmp.lty <- c(ghyp.lty, "dashed")
                legend(location, legend = tmp.text, col = tmp.col,
                       lty = tmp.lty, cex = legend.cex)
            }else{
                legend(location, legend = ghyp.name(x, abbr = TRUE, skew.attr = TRUE), col = ghyp.col,
                       lty = ghyp.lty, cex = legend.cex)
            }
        }
    }
}
### <---------------------------------------------------------------------->
setMethod("hist", signature(x = "ghyp"), hist.ghyp)
### <---------------------------------------------------------------------->


### <======================================================================>
"lines.ghyp" <- function(x, range = qghyp(c(0.001, 0.999), x),
                           length = 1000, ...)
{
    .test.ghyp(x, "univariate")
    if(length(range) > 2){
        x.seq <- range
    }else{
        x.seq <- seq(min(range), max(range), length = length)
    }
    ghyp.density <- dghyp(x.seq, x)
    lines(x.seq, ghyp.density, ...)
}
### <---------------------------------------------------------------------->
setMethod("lines", signature(x = "ghyp"), lines.ghyp)
### <---------------------------------------------------------------------->


### <======================================================================>
"logLik.mle.ghyp" <- function(object, ...)
{
    ghyp.models <- list(...)
    ## Internal function; test whether all objects are of type mle.ghyp
    test.class.mle.ghyp <- function(x){
        if(!is(x, "mle.ghyp")){
            stop("Some of the objects are not of class 'mle.ghyp'!")
        }
    }

    sapply(ghyp.models, test.class.mle.ghyp)

    ghyp.models <- c(object, ghyp.models)
    tmp.logLik <- sapply(ghyp.models, function(z)ghyp.fit.info(z)$logLikelihood)
    return(tmp.logLik)
}
### <---------------------------------------------------------------------->
setMethod("logLik", signature(object = "mle.ghyp"), logLik.mle.ghyp)
### <---------------------------------------------------------------------->



### <======================================================================>
"mean.ghyp" <- function(x)
{
    return(x@expected.value)
}
### <---------------------------------------------------------------------->
setMethod("mean", signature(x = "ghyp"), mean.ghyp)
### <---------------------------------------------------------------------->



### <======================================================================>
"pairs.ghyp" <- function(x, data = ghyp.data(x), main = "'ghyp' pairwise plot",
                           nbins = 30, qq = TRUE, gaussian = TRUE,
                           hist.col = c("white", topo.colors(40)),
                           spline.points = 150, root.tol = .Machine$double.eps^0.5,
                           rel.tol = root.tol, abs.tol = root.tol^1.5, ...)
{
    .test.ghyp(x, case = "multivariate")

    data <- .check.data(data, case = "mv", na.rm = T, fit = TRUE, dim = x@dimension)

    nc <- length(x@mu)

    ## local function for the purpose to plot axes (copied from pairs.default)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main,oma, ...) {
        if (side%%2 == 1){
            Axis(x, side = side, xpd = NA, ...)
        }else{
            Axis(y, side = side, xpd = NA, ...)
        }
    }

    ## If not defined: define a sensible title
    if(main=="'ghypmv' pairwise plot."){
        if(qq & !is.null(colnames(data))){
            main <- paste(main, "\nColnames: ",
                          paste(colnames(data), collapse = ", "), sep="")
        }
    }

    tmp.par <- par()
    plot.new()

    par(mfrow = dim(x@sigma), mar = c(0, 0, 0, 0) + .1, oma =   c(0, 0, 4, 0) + 3)

    on.exit(suppressWarnings(par(tmp.par)))

    for(i in 1:nc){
        for(j in 1:nc){
            par(mfg=c(i,j))
            if(i==j){
                if(qq){
                    qqghyp(x[i], data[,i], main = "", plot.legend = gaussian,
                           gaussian = gaussian, xaxt = "n", yaxt = "n",
                           spline.points = spline.points, root.tol = root.tol,
                           rel.tol = rel.tol, abs.tol = abs.tol, ...)

                }else{
                    plot(NA, xlim = c(0,1), ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "",...)
                    text(0.5, 0.5, labels = colnames(data)[i])
                }
            }else{

                if(i > j){
                    ## 2 dimensional histogramm from package 'gplots'
                    hist2d(data[, c(j, i)], nbins = nbins, col = hist.col, xaxt = "n", yaxt = "n",
                           xlab = "", ylab = "", axes = T, frame.plot = T, bty = "o")
                    box()
                }else{
                    plot((data[, c(j, i)]), xaxt = "n", yaxt = "n", xlab = "", ylab = "", ...)
                }
            }

            if (j == 1 && !(i%%2))
                Axis(data[, i], side=2,...)

            if (i == 1 && !(j%%2))
                Axis(data[, j], side=3, ...)

            if (i == nc && (j%%2))
                Axis(data[, j], side=1, ...)
            if (j == nc && (i%%2))
                Axis(data[, i], side=4, ...)
        }
    }
    title(main = main, outer = T)
}
### <---------------------------------------------------------------------->
setMethod("pairs", signature(x = "ghyp"), pairs.ghyp)
### <---------------------------------------------------------------------->



### <======================================================================>
"plot.ghyp" <- function(x, range = qghyp(c(0.001, 0.999), x),
                          length = 1000, ...)
{
    .test.ghyp(x, "univariate")
    if(length(range) > 2){
        x.seq <- range
    }else{
        x.seq <- seq(min(range), max(range), length = length)
    }
    ghyp.density <- dghyp(x.seq, x)
    plot(x.seq, ghyp.density, ...)
}
### <---------------------------------------------------------------------->
setMethod("plot", signature(x = "ghyp",y = "missing"), plot.ghyp)
### <---------------------------------------------------------------------->



### <======================================================================>
scale.ghyp <- function(x, center = TRUE, scale = TRUE)
{
    if(scale & any(!is.finite(vcov(x))))
    {
        stop("Object must have finite variance if 'scale == TRUE'!")
    }

    if(center){
        x <- transform(x, summand = - mean(x))
    }

    if(scale){
        if(.is.univariate(x)){
            x <- transform(x, multiplier = 1/sqrt(vcov(x)))
        }else{
            multiplier <- chol(solve(vcov(x)))
            x <- transform(x, multiplier = multiplier)
        }

    }
    return(x)
}
### <---------------------------------------------------------------------->
setMethod("scale", signature(x = "ghyp"), scale.ghyp)
### <---------------------------------------------------------------------->



### <======================================================================>
"show.ghyp" <- function(object)
{
    cat(ghyp.name(object, abbr = FALSE, skew.attr = TRUE), "Distribution:\n", sep = " ")
    cat("\nParameters:\n")
    if(ghyp.dim(object) > 1){           # Multivariate case
        param <- coef(object)
        if(object@parametrization == "alpha.delta"){
            if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "t"){
                ## Student-t  ->  alpha^2 == beta' Delta beta
                param.uv <- unlist(param[c(1, 3)])
            }else if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "VG"){
                ## VG  ->  delta == 0
                param.uv <- unlist(param[1:2])
            }else if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "ghyp"){
                param.uv <- unlist(param[1:3])
            }else{
                ## hyp or NIG -> lambda is constant
                param.uv <- unlist(param[2:3])
            }
            print(param.uv)
            cat("\nmu:\n")
            print(param$mu)
            cat("\nDelta:\n")
            print(param$Delta)
            cat("\nbeta:\n")
            print(param$beta)
        }else if(.is.gaussian(object)){
            cat("\nmu:\n")
            print(param$mu)
            cat("\nsigma:\n")
            print(param$sigma)
        }else{
            if(object@parametrization == "chi.psi"){
                param.uv <- unlist(param[1:3])
            }else{
                param.uv <- unlist(param[1:2])
            }
            if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "t"){
                ## Student-t  ->  alpha.bar == 0
                if(object@parametrization == "chi.psi"){
                    param.uv <- param.uv[c("lambda", "chi")]
                    param.uv <- c(nu = unname(-2 * param.uv["lambda"]), param.uv["chi"])
                }else{
                    param.uv <- param.uv["nu"]
                }
            }else if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "VG"){
                ## VG  ->  alpha.bar == 0 or chi == 0
                if(object@parametrization == "chi.psi"){
                    param.uv <- param.uv[c(1, 3)]
                }else if(object@parametrization == "alpha.bar"){
                    param.uv <- param.uv[1]
                }
            }else if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "ghyp"){
                ##param.uv <- param.uv
            }else{
                ## hyp or NIG -> lambda is constant
                param.uv <- param.uv[-1]
            }
            print(param.uv)
            cat("\nmu:\n")
            print(param$mu)
            cat("\nsigma:\n")
            print(param$sigma)
            cat("\ngamma:\n")
            print(param$gamma)
        }
    }else{                              # Univariate case
        param <- unlist(coef(object))
        ## alpha.delta-parametrization
        if(.is.gaussian(object)){
            param <- c(param["mu"], param["sigma"])
        }else if(object@parametrization == "alpha.delta"){
            if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "t"){
                ## Student-t  ->  alpha^2 == beta^2
                param <- param[-2]
            }else if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "VG"){
                ## VG  ->  delta == 0
                param <- param[-5]
            }else if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "ghyp"){
                                        # param <- param
            }else{
                                        # hyp or NIG -> lambda is constant
                param <- param[-1]
            }
        }else {
            ## chi.psi or alpha.bar-parametrization
            if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "t"){
                ## Student-t  ->  alpha.bar == 0
                if(object@parametrization == "chi.psi"){
                    param <- c(nu = -2 * unname(param[1]), param[-c(1, 3)])
                }else{
                    param <- param[-c(1, 3)]
                }
            }else if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "VG"){
                ## VG  ->  alpha.bar == 0 or chi == 0
                param <- param[-2]
            }else if(ghyp.name(object, abbr = TRUE, skew.attr = FALSE) == "ghyp"){
                ## param <- param
            }else{
                ## hyp or NIG -> lambda is constant
                param <- param[-1]
            }
        }
        print(param)
    }

}
### <---------------------------------------------------------------------->
setMethod("show", signature(object = "ghyp"), show.ghyp)
### <---------------------------------------------------------------------->



### <======================================================================>
"show.mle.ghyp" <- function(object)
{
    if(!object@converged && !is.na(object@n.iter)){
        cat("Warning: fitting procedure did not converge!\n\n")
    }else if(!object@converged && is.na(object@n.iter)){
        cat("Error: fitting procedure crashed! Use 'summary' to see the message!\n\n")
    }

    callNextMethod()

    cat("\nlog-likelihood:\n", logLik(object), "\n\n", sep = "")

    cat("\nCall:\n", deparse(object@call), "\n\n", sep = "")
}
### <---------------------------------------------------------------------->
setMethod("show", signature(object = "mle.ghyp"), show.mle.ghyp)
### <---------------------------------------------------------------------->


### <======================================================================>
"summary.mle.ghyp" <- function(object)
{
    if(!object@converged && !is.na(object@n.iter)){
        cat("Warning: fitting procedure did not converge!\n\n")
    }else if(!object@converged && is.na(object@n.iter)){
        cat("Error: fitting procedure crashed! See error message below!\n\n")
    }

    show.ghyp(object)

    cat("\nCall:\n", deparse(object@call), "\n\n", sep = "")

    cat("Optimization information:\n")
    cat("log-Likelihood:               ", logLik(object), "\n")
    cat("AIC:                          ", AIC(object),    "\n")

    object.dim <- ghyp.dim(object)
    if(.is.gaussian(object)){
        nbr.fitted.params <- unname(object.dim + object.dim/2 * (object.dim + 1) *
                                    object@fitted.params["sigma"])

        names.fitted.param <- paste(names(object@fitted.params[object@fitted.params]),collapse=", ")
    }else{
        nbr.fitted.params <- unname(sum(object@fitted.params[c("alpha.bar", "lambda")]) +
                                    object.dim * sum(object@fitted.params[c("mu", "gamma")]) +
                                    object.dim/2 * (object.dim + 1) * object@fitted.params["sigma"])

        names.fitted.param <- paste(names(object@fitted.params[object@fitted.params]),collapse=", ")
    }
    cat("Fitted parameters:             " , names.fitted.param,   ";  (Number: ", nbr.fitted.params, ")\n", sep = "")
    cat("Number of iterations:         "  , object@n.iter,        "\n")
    cat("Converged:                    "  , object@converged,     "\n")
    if(!object@converged){
        cat("Error code:                   ", object@error.code,    "\n")
        cat("Error message:                ", object@error.message, "\n")
    }
}
### <---------------------------------------------------------------------->
setMethod("summary", signature(object = "mle.ghyp"), summary.mle.ghyp)
### <---------------------------------------------------------------------->



### <======================================================================>
"transform.ghyp" <- function(`_data`, summand, multiplier)
{
    ##  show(`_data`)
    ##  print(substitute(`_data`))
    ##  return(TRUE)
    object <- `_data`
    .test.ghyp(object, case = "ghyp")

    if(!missing(multiplier))
    {
        if(any(!is.finite(multiplier)))
        {
            stop("All elements of 'multiplier' must be finite!")
        }
    }

    if(!missing(summand))
    {
        if(any(!is.finite(summand)))
        {
            stop("All elements of 'summand' must be finite!")
        }
    }

    if(!missing(multiplier)){
        tmp.multiplier <- as.matrix(multiplier)
        if(max(dim(tmp.multiplier)) > ghyp.dim(object))
            stop("Extending the dimension of a ghyp object is not allowed!",
                 "That is, ncol(multiplier) and length(summand) must be <= dimension!")

    }

    ## There are 4 cases:
    ## (1) 'summand' is missing: Set 'summand' to 0
    ## (2) 'multiplier' is missing: Set 'multiplier' to the unit matrix
    ## (3) 'multiplier' and 'summand' are passed:
    ##     The dimension of 'multiplier' is k x n : 'n' must be the dimension
    ##                                              ghyp object and 'k' must be
    ##                                              the length of 'summand'
    ## (4) 'multiplier' and 'summand' are missing: -> stop

    ## When the map is from R^n to R (i.e. multivariate to univariate), sigma acutally
    ## changes its interpretation from variance to standard deviation

    if(missing(summand) & missing(multiplier)){
        ## case (3)
        stop("No arguments submitted!")
    }else if(missing(summand)){
        ## case (1)
        if(!is.matrix(multiplier)){
            multiplier <- matrix(multiplier, nrow = 1)
        }
        if(ncol(multiplier) != ghyp.dim(object)){
            stop("Dimension of multiplier must be ",
                 "n x ", ghyp.dim(object), "!")
        }
        summand <- rep(0, nrow(multiplier))
    }else if(missing(multiplier)){
        ## case (2)
        summand <- as.vector(summand)
        if(length(summand) != ghyp.dim(object)){
            stop("Dimension of summand must be ", ghyp.dim(object), "!")
        }
        multiplier <- diag(rep(1, ghyp.dim(object)))
    }else{
        ## case (4)
        summand <- as.vector(summand)
        if(!is.matrix(multiplier)){
            multiplier <- matrix(multiplier, nrow = 1)
        }
        if(ncol(multiplier) != ghyp.dim(object)){
            stop("Dimension of multiplier must be ",
                 "n x ", ghyp.dim(object), "!")
        }
        if(length(summand) != nrow(multiplier)){
            stop("Dimension mimatch: length(summand) must be equal to nrow(multiplier)!")
        }
    }

    if(any(diag(multiplier) == 0)){
        stop("Diagonal elements of multiplier must not be '0'!")
    }

    sigma <- multiplier %*% object@sigma %*% t(multiplier)
    if(ncol(sigma) == 1){
        sigma <- as.vector(sigma)
        if(.is.univariate(object)){
            sigma <- sigma * object@sigma # = multiplier^2 * sigma^2 -> univariate case
        }
    }

    if(length(summand) == 1){
        sigma <- sqrt(sigma)
    }

    if(is.null(ghyp.data(object))){
        transformed.data <- NULL
    }else{
        transformed.data <- t(multiplier %*% t(ghyp.data(object)) +
            t(matrix(rep(summand, nrow(as.matrix(ghyp.data(object)))),
                     byrow = TRUE, ncol = length(summand))))
    }

    ## Return the transformed object in the same parametrization
    if(object@parametrization == "alpha.bar"){
        if(.is.student.t(object)){
            return(student.t(nu = -2 * object@lambda,
                             mu = as.vector(multiplier %*% object@mu + summand),
                             sigma = sigma,
                             gamma = as.vector(multiplier %*% object@gamma),
                             data = transformed.data))
        }else{
            return(ghyp(lambda = object@lambda, alpha.bar = object@alpha.bar,
                        mu = as.vector(multiplier %*% object@mu + summand),
                        sigma = sigma, gamma = as.vector(multiplier %*% object@gamma),
                        data = transformed.data))
        }
    }else if(object@parametrization == "chi.psi"){
        return(ghyp(lambda = object@lambda, chi = object@chi, psi = object@psi,
                    mu = as.vector(multiplier %*% object@mu + summand),
                    sigma = sigma, gamma = as.vector(multiplier %*% object@gamma),
                    data = transformed.data))
    }else if(object@parametrization == "alpha.delta"){
        tmp.obj <- ghyp(lambda = object@lambda, chi = object@chi, psi = object@psi,
                        mu = as.vector(multiplier %*% object@mu + summand),
                        sigma = sigma, gamma = as.vector(multiplier %*% object@gamma))
        return(do.call("ghyp.ad", c(coef(tmp.obj, type = "alpha.delta"), list(data = transformed.data))))
    }else if(object@parametrization == "Gaussian"){
        return(gauss(mu = as.vector(multiplier %*% object@mu + summand),
                     sigma = sigma, data = transformed.data))
    }else{
        stop("Unknown parametrization: ", object@parametrization)
    }
}
### <---------------------------------------------------------------------->
setMethod("transform", signature(`_data` = "ghyp"), transform.ghyp)
### <---------------------------------------------------------------------->


### <======================================================================>
"vcov.ghyp" <- function(object)
{
    if(.is.univariate(object)){
        return(as.vector(object@variance))
    }else{
        return(object@variance)
    }
}
### <---------------------------------------------------------------------->
setMethod("vcov", signature(object = "ghyp"), vcov.ghyp)
### <---------------------------------------------------------------------->
