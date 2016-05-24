fitspatgev <- function(data, covariables, loc.form, scale.form, shape.form,
                       temp.cov = NULL, temp.form.loc = NULL, temp.form.scale = NULL,
                       temp.form.shape = NULL, ..., start, control = list(maxit = 10000),
                       method = "Nelder", warn = TRUE, corr = FALSE){

    n.site <- ncol(data)
    n.obs <- nrow(data)

    if (n.site != nrow(covariables))
        stop("'data' and 'covariates' doesn't match")

    use.temp.cov <- c(!is.null(temp.form.loc), !is.null(temp.form.scale), !is.null(temp.form.shape))

    if (any(use.temp.cov) && (n.obs != nrow(temp.cov)))
        stop("'data' and 'temp.cov' doesn't match")

    if (any(use.temp.cov) && is.null(temp.cov))
        stop("'temp.cov' must be supplied if at least one temporal formula is given")

    ##With our notation, formula must be of the form y ~ xxxx
    loc.form <- update(loc.form, y ~ .)
    scale.form <- update(scale.form, y ~ .)
    shape.form <- update(shape.form, y ~ .)

    if (use.temp.cov[1])
        temp.form.loc <- update(temp.form.loc, y ~. + 0)

    if (use.temp.cov[2])
        temp.form.scale <- update(temp.form.scale, y ~. + 0)

    if (use.temp.cov[3])
        temp.form.shape <- update(temp.form.shape, y ~. + 0)

    loc.model <- modeldef(covariables, loc.form)
    scale.model <- modeldef(covariables, scale.form)
    shape.model <- modeldef(covariables, shape.form)

    loc.dsgn.mat <- loc.model$dsgn.mat
    scale.dsgn.mat <- scale.model$dsgn.mat
    shape.dsgn.mat <- shape.model$dsgn.mat

    loc.pen.mat <- loc.model$pen.mat
    scale.pen.mat <- scale.model$pen.mat
    shape.pen.mat <- shape.model$pen.mat

    loc.penalty <- loc.model$penalty.tot
    scale.penalty <- scale.model$penalty.tot
    shape.penalty <- shape.model$penalty.tot

    ##The total number of parameters to be estimated for each GEV
    ##parameter
    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)

    ##The number of ``purely parametric'' parameters to estimate i.e. we
    ##do not consider the weigths given to each basis function
    n.pparloc <- loc.model$n.ppar
    n.pparscale <- scale.model$n.ppar
    n.pparshape <- shape.model$n.ppar

    loc.names <- paste("locCoeff", 1:n.loccoeff, sep="")
    scale.names <- paste("scaleCoeff", 1:n.scalecoeff, sep="")
    shape.names <- paste("shapeCoeff", 1:n.shapecoeff, sep="")

    ##Do the same for the temporal regression coefficients
    if (use.temp.cov[1]){
        temp.model.loc <- modeldef(temp.cov, temp.form.loc)
        temp.dsgn.mat.loc <- temp.model.loc$dsgn.mat
        temp.pen.mat.loc <- temp.model.loc$pen.mat
        temp.penalty.loc <- temp.model.loc$penalty.tot
        n.tempcoeff.loc <- ncol(temp.dsgn.mat.loc)
        n.ppartemp.loc <- temp.model.loc$n.ppar
        temp.names.loc <- paste("tempCoeffLoc", 1:n.tempcoeff.loc, sep="")
    }

    else {
        temp.model.loc <- temp.dsgn.mat.loc <- temp.pen.mat.loc <- temp.names.loc <- NULL
        n.tempcoeff.loc <- n.ppartemp.loc <- temp.penalty.loc <- 0
    }

    if (use.temp.cov[2]){
        temp.model.scale <- modeldef(temp.cov, temp.form.scale)
        temp.dsgn.mat.scale <- temp.model.scale$dsgn.mat
        temp.pen.mat.scale <- temp.model.scale$pen.mat
        temp.penalty.scale <- temp.model.scale$penalty.tot
        n.tempcoeff.scale <- ncol(temp.dsgn.mat.scale)
        n.ppartemp.scale <- temp.model.scale$n.ppar
        temp.names.scale <- paste("tempCoeffScale", 1:n.tempcoeff.scale, sep="")
    }

    else {
        temp.model.scale <- temp.dsgn.mat.scale <- temp.pen.mat.scale <- temp.names.scale <- NULL
        n.tempcoeff.scale <- n.ppartemp.scale <- temp.penalty.scale <- 0
    }

    if (use.temp.cov[3]){
        temp.model.shape <- modeldef(temp.cov, temp.form.shape)
        temp.dsgn.mat.shape <- temp.model.shape$dsgn.mat
        temp.pen.mat.shape <- temp.model.shape$pen.mat
        temp.penalty.shape <- temp.model.shape$penalty.tot
        n.tempcoeff.shape <- ncol(temp.dsgn.mat.shape)
        n.ppartemp.shape <- temp.model.shape$n.ppar
        temp.names.shape <- paste("tempCoeffShape", 1:n.tempcoeff.shape, sep="")
    }

    else {
        temp.model.shape <- temp.dsgn.mat.shape <- temp.pen.mat.shape <- temp.names.shape <- NULL
        n.tempcoeff.shape <- n.ppartemp.shape <- temp.penalty.shape <- 0
    }

    param <- c(loc.names, scale.names, shape.names, temp.names.loc, temp.names.scale,
               temp.names.shape)

    nllik <- function(x) x

    body(nllik) <- parse(text = paste("-.C('spatgevlik', as.double(data), as.double(covariables),
 as.integer(n.site), as.integer(n.obs), as.double(loc.dsgn.mat), as.double(loc.pen.mat),
 as.integer(n.loccoeff), as.integer(n.pparloc), as.double(loc.penalty),
 as.double(scale.dsgn.mat), as.double(scale.pen.mat), as.integer(n.scalecoeff),
 as.integer(n.pparscale), as.double(scale.penalty), as.double(shape.dsgn.mat),
 as.double(shape.pen.mat), as.integer(n.shapecoeff), as.integer(n.pparshape),
 as.double(shape.penalty), as.integer(use.temp.cov), as.double(temp.dsgn.mat.loc),
 as.double(temp.pen.mat.loc), as.integer(n.tempcoeff.loc), as.integer(n.ppartemp.loc),
 as.double(temp.penalty.loc), as.double(temp.dsgn.mat.scale), as.double(temp.pen.mat.scale),
 as.integer(n.tempcoeff.scale), as.integer(n.ppartemp.scale), as.double(temp.penalty.scale),
 as.double(temp.dsgn.mat.shape), as.double(temp.pen.mat.shape), as.integer(n.tempcoeff.shape),
 as.integer(n.ppartemp.shape), as.double(temp.penalty.shape),",
                         paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                         paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                         paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                         paste("as.double(c(", paste(temp.names.loc, collapse = ","), ")), "),
                         paste("as.double(c(", paste(temp.names.scale, collapse = ","), ")), "),
                         paste("as.double(c(", paste(temp.names.shape, collapse = ","), ")), "),
                         "dns = double(1), PACKAGE = 'SpatialExtremes', NAOK = TRUE)$dns"))


    ##Define the formal arguments of the function
    form.nllik <- NULL
    for (i in 1:length(param))
        form.nllik <- c(form.nllik, alist(a=))

    names(form.nllik) <- param
    formals(nllik) <- form.nllik

    if (missing(start)){
        loc <- scale <- shape <- rep(0, n.site)

        for (i in 1:n.site){
            gev.param <- gevmle(data[,i])
            loc[i] <- gev.param["loc"]
            scale[i] <- gev.param["scale"]
            shape[i] <- gev.param["shape"]
        }

        locCoeff <- loc.model$init.fun(loc)
        scaleCoeff <- scale.model$init.fun(scale)
        shapeCoeff <- shape.model$init.fun(shape)

        locCoeff[is.na(locCoeff)] <- 0
        scaleCoeff[is.na(scaleCoeff)] <- 0
        shapeCoeff[is.na(shapeCoeff)] <- 0

        ##To be sure that the scale parameter is always positive at starting
        ##values
        scales.hat <- scale.model$dsgn.mat %*% scaleCoeff

        if (any(scales.hat <= 0))
            scaleCoeff[1] <- scaleCoeff[1] - 1.001 * min(scales.hat)

        names(locCoeff) <- loc.names
        names(scaleCoeff) <- scale.names
        names(shapeCoeff) <- shape.names

        if (use.temp.cov[1]){
            tempCoeff.loc <- rep(0, n.tempcoeff.loc)
            names(tempCoeff.loc) <- temp.names.loc
        }

        else
            tempCoeff.loc <- NULL

        if (use.temp.cov[2]){
            tempCoeff.scale <- rep(0, n.tempcoeff.scale)
            names(tempCoeff.scale) <- temp.names.scale
        }

        else
            tempCoeff.scale <- NULL

        if (use.temp.cov[3]){
            tempCoeff.shape <- rep(0, n.tempcoeff.shape)
            names(tempCoeff.shape) <- temp.names.shape
        }

        else
            tempCoeff.shape <- NULL

        start <- as.list(c(locCoeff, scaleCoeff, shapeCoeff, tempCoeff.loc,
                           tempCoeff.scale, tempCoeff.shape))
        start <- start[!(param %in% names(list(...)))]

    }

    if (!length(start))
        stop("there are no parameters left to maximize over")

    nm <- names(start)
    l <- length(nm)
    f <- formals(nllik)
    names(f) <- param
    m <- match(nm, param)

    if(any(is.na(m)))
        stop("'start' specifies unknown arguments")

    formals(nllik) <- c(f[m], f[-m])
    nllh <- function(p, ...) nllik(p, ...)

    if(l > 1)
        body(nllh) <- parse(text = paste("nllik(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))

    fixed.param <- list(...)[names(list(...)) %in% param]

    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")

    start.arg <- c(list(p = unlist(start)), fixed.param)

    init.lik <- do.call("nllh", start.arg)
    if (warn && (init.lik >= 1.0e6))
        warning("negative log-likelihood is infinite at starting values")

    if (method == "nlminb"){
        start <- as.numeric(start)
        opt <- nlminb(start, nllh, ..., control = control)
        opt$counts <- opt$evaluations
        opt$value <- opt$objective
        names(opt$par) <- nm
    }

    else if (method == "nlm"){
        start <- as.numeric(start)
        opt <- nlm(nllh, start, ...)
        opt$counts <- opt$iterations
        names(opt$counts) <- "function"
        opt$value <- opt$minimum
        opt$par <- opt$estimate
        names(opt$par) <- nm

        if (opt$code <= 2){
            opt$convergence <- 0
            opt$message <- NULL
        }

        if (opt$code > 2){
            opt$convergence <- 1
            opt$message <- paste("nlm error code", opt$code)
        }
    }

    else
        opt <- optim(start, nllh, ..., method = method, control = control)

    if ((opt$convergence != 0) || (opt$value >= 1.0e6)){
        if (warn)
            warning("optimization may not have succeeded")
    }

    else
        opt$convergence <- "successful"

    param.names <- param
    param <- c(opt$par, unlist(fixed.param))
    param <- param[param.names]

    std.err <- .spatgevstderr(param, data, loc.dsgn.mat,scale.dsgn.mat, shape.dsgn.mat,
                              temp.dsgn.mat.loc, temp.dsgn.mat.scale, temp.dsgn.mat.shape,
                              use.temp.cov, fixed.param = names(fixed.param),
                              param.names = param.names)

    opt$hessian <- std.err$hessian
    var.score <- std.err$var.score
    ihessian <- try(solve(opt$hessian), silent = TRUE)

    if(!is.matrix(ihessian) || any(is.na(var.score))){
        if (warn)
            warning("Observed information matrix is singular. No standard error will be computed.")

        std.err.type <- "none"
    }

    else{
        std.err.type <- "yes"
        var.cov <- ihessian %*% var.score %*% ihessian
        std.err <- diag(var.cov)

        std.idx <- which(std.err <= 0)
        if(length(std.idx) > 0){
            if (warn)
                warning("Some (observed) standard errors are negative;\n passing them to NA")

            std.err[std.idx] <- NA
        }


        std.err <- sqrt(std.err)

        if(corr) {
            .mat <- diag(1/std.err, nrow = length(std.err))
            corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
            diag(corr.mat) <- rep(1, length(std.err))
        }

        else
            corr.mat <- NULL

        colnames(var.cov) <- colnames(ihessian) <- rownames(var.cov) <-
            rownames(ihessian) <- names(std.err) <- nm
    }

    if (std.err.type == "none")
        std.err <- std.err.type <- corr.mat <- var.cov <- ihessian <-
            var.score <- NULL

    ans <- list(fitted.values = opt$par, param = param, std.err = std.err, var.cov = var.cov,
                counts = opt$counts, message = opt$message, covariables = covariables,
                logLik = -opt$value, loc.form = loc.form, scale.form = scale.form,
                shape.form = shape.form, convergence = opt$convergence, nllh = nllh,
                deviance = 2 * opt$value, ihessian = ihessian, var.score = var.score,
                data = data, fixed = unlist(fixed.param), hessian = opt$hessian,
                use.temp.cov = use.temp.cov)

    class(ans) <- "spatgev"
    return(ans)
}
