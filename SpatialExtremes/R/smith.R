##This file contains all the function to fit the max-stable
##characterisation of Smith


##This functions fits the model without any spatial structure for the
##GEV parameters. Thus, each GEV parameters are estimated at each
##location. However, if fit.marge = FALSE, observation are supposed to
##be unit Frechet and only the covariance matrix is estimated.
smithfull <- function(data, coord, start, fit.marge = FALSE, iso = TRUE,
                      ..., warn = TRUE, method = "BFGS",
                      control = list(),
                      corr = FALSE, weights = NULL, check.grad = FALSE){
    ##data is a matrix with each column corresponds to one location
    ##coord is a matrix giving the coordinates (1 row = 1 station)
    n.site <- ncol(data)
    n.obs <- nrow(data)
    dist.dim <- ncol(coord)
    n.pair <- n.site * (n.site - 1) / 2

    distVec <- distance(coord, vec = TRUE)
    weighted <- !is.null(weights)

    if (!weighted)
        ##Set the weights to 0 as they won't be used anyway
        weights <- 0

    if (iso)
        param <- "cov"

    else{
        if (dist.dim == 2)
            param <- c("cov11", "cov12", "cov22")

        else if (dist.dim == 3)
            param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33")
    }

    if (fit.marge){
        loc.names <- paste("loc", 1:n.site, sep="")
        scale.names <- paste("scale", 1:n.site, sep="")
        shape.names <- paste("shape", 1:n.site, sep="")
        param <- c(param, loc.names, scale.names, shape.names)
    }

    else
        loc.names <- scale.names <- shape.names <- rep(1, n.site)

    ##First create a "void" function
    nplk <- function(x) x

    ##And define the "body" of the function as the number of parameters
    ##to estimate depends on n.site
    if (dist.dim == 2){
        if (iso)
            body(nplk) <- parse(text = paste("-.C('smithfull', as.double(data), as.double(distVec), as.integer(n.site), as.integer(n.obs), as.integer(weighted), as.double(weights),",
                                paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                                "as.double(cov), as.double(0.0), as.double(cov), fit.marge, dns = double(1), PACKAGE = 'SpatialExtremes', NAOK = TRUE)$dns"))

        else
            body(nplk) <- parse(text = paste("-.C('smithfull', as.double(data), as.double(distVec), as.integer(n.site), as.integer(n.obs), as.integer(weighted), as.double(weights),",
                                paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                                "as.double(cov11), as.double(cov12), as.double(cov22), fit.marge, dns = double(1), PACKAGE = 'SpatialExtremes', NAOK = TRUE)$dns"))
    }

    if (dist.dim == 3){
        if (iso)
            body(nplk) <- parse(text = paste("-.C('smithfull3d', as.double(data), as.double(distVec), as.integer(n.site), as.integer(n.obs), as.integer(weighted), as.double(weights),",
                                paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                                "as.double(cov), as.double(0.0), as.double(0.0), as.double(cov), as.double(0.0), as.double(cov), fit.marge, dns = double(1), PACKAGE = 'SpatialExtremes, NAOK = TRUE')$dns"))

        else
            body(nplk) <- parse(text = paste("-.C('smithfull3d', as.double(data), as.double(distVec), as.integer(n.site), as.integer(n.obs), as.integer(weighted), as.double(weights),",
                                paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                                "as.double(cov11), as.double(cov12), as.double(cov13), as.double(cov22), as.double(cov23), as.double(cov33), fit.marge, dns = double(1), PACKAGE = 'SpatialExtremes, NAOK = TRUE')$dns"))
    }

    fixed.param <- list(...)[names(list(...)) %in% param]

    ##Define the formal arguments of the function
    form.nplk <- NULL
    for (i in 1:length(param))
        form.nplk <- c(form.nplk, alist(a=))

    names(form.nplk) <- param
    formals(nplk) <- form.nplk

    if (missing(start)) {

        start <- list()
        if (fit.marge){
            locs <- scales <- rep(NA, n.site)
            shapes <- rep(0, n.site)

            for (i in 1:n.site){
                marg.param <- gevmle(data[,i])
                locs[i] <- marg.param["loc"]
                scales[i] <- marg.param["scale"]
                shapes[i] <- marg.param["shape"]
            }

            start <- as.list(unlist(list(loc = locs, scale = scales, shape = shapes)))
        }

        if (length(fixed.param) > 0){
            args <- c(list(data = data, coord = coord, marge = "emp", iso = iso, fixed.param))
            cov.start <- do.call("fitcovmat", args)$param
        }

        else
            cov.start <- fitcovmat(data, coord, marge = "emp", iso = iso)$param

        if (iso){
            cov.start <- cov.start[1]
            names(cov.start) <- "cov"
        }

        start <- c(as.list(cov.start), start)
        start <- start[!(param %in% names(list(...)))]
    }


    if (!is.list(start))
        stop("'start' must be a named list")

    if (!length(start))
        stop("there are no parameters left to maximize over")

    nm <- names(start)
    l <- length(nm)
    f <- formals(nplk)
    names(f) <- param
    m <- match(nm, param)

    if(any(is.na(m)))
        stop("'start' specifies unknown arguments")

    formals(nplk) <- c(f[m], f[-m])
    nllh <- function(p, ...) nplk(p, ...)

    if(l > 1)
        body(nllh) <- parse(text = paste("nplk(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))

    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")

    start.arg <- c(list(p = unlist(start)), fixed.param)

    init.lik <- do.call("nllh", start.arg)
    if (warn && (init.lik >= 1.0e15))
        warning("negative log-likelihood is infinite at starting values")

    if (method == "nlminb"){
        start <- as.numeric(start)
        opt <- nlminb(start, nllh, ..., control = control)
        opt$counts <- opt$evaluations
        opt$value <- opt$objective
        names(opt$par) <- nm

        if ((opt$convergence != 0) || (opt$value >= 1.0e15)) {
            if (warn)
                warning("optimization may not have succeeded")
        }

        if (opt$convergence == 0)
            opt$convergence <- "successful"
    }

    if (method == "nlm"){
        start <- as.numeric(start)
        opt <- nlm(nllh, start, ...)
        opt$counts <- opt$iterations
        names(opt$counts) <- "function"
        opt$value <- opt$minimum
        opt$par <- opt$estimate
        names(opt$par) <- nm

        if (opt$code <= 2)
            opt$convergence <- "sucessful"

        if (opt$code == 3)
            opt$convergence <- "local minimum or 'steptol' is too small"

        if (opt$code == 4)
            opt$convergence <- "iteration limit reached"

        if (opt$code == 5)
            opt$convergence <- "optimization failed"
    }

    if (!(method %in% c("nlm", "nlminb"))){
        opt <- optim(start, nllh, ..., method = method, control = control)

        if ((opt$convergence != 0) || (opt$value >= 1.0e15)) {
            if (warn)
                warning("optimization may not have succeeded")

            if (opt$convergence == 1)
                opt$convergence <- "iteration limit reached"
        }

        else opt$convergence <- "successful"
    }

    if (opt$value == init.lik){
        if (warn)
            warning("optimization stayed at the starting values.")

        opt$convergence <- "Stayed at start. val."
    }

    param.names <- param
    param <- c(opt$par, unlist(fixed.param))
    param <- param[param.names]

    if (iso){
        if (dist.dim == 2){
            param <- c(param["cov"], 0, param["cov"], param[-1])
            names(param)[1:3] <- c("cov11", "cov12", "cov22")
        }

        else{
            param <- c(param["cov"], 0, 0, param["cov"], 0, param["cov"],
                       param[-1])
            names(param)[1:6] <- c("cov11", "cov12", "cov13", "cov22", "cov23",
                                   "cov33")
        }
    }

    ##Reset the weights to their original values
    if ((length(weights) == 1) && (weights == 0))
        weights <- NULL


    std.err <- .smithstderr(param, data, distVec, as.double(0), as.double(0),
                            as.double(0), as.double(0), as.double(0), as.double(0),
                            rep(FALSE, 3), fit.marge = fit.marge,
                            fixed.param = names(fixed.param),
                            param.names = param.names, iso = iso, weights = weights)

    if (check.grad)
        print(round(rbind(numerical = -opt$grad, analytical = std.err$grad), 3))

    opt$hessian <- std.err$hessian
    var.score <- std.err$var.score
    ihessian <- try(solve(opt$hessian), silent = TRUE)

    if(!is.matrix(ihessian) || any(is.na(var.score))){
        if (warn)
            warning("Observed information matrix is singular. No standard error will be computed.")

        std.err.type <- "none"
    }

    else {
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


    if (std.err.type == "none"){
        std.err <- std.err.type <- corr.mat <- NULL
        var.cov <- ihessian <- var.score <- NULL
    }

    if (dist.dim == 2)
        Sigma <- matrix(c(param["cov11"], param["cov12"], param["cov12"],
                          param["cov22"]), 2, 2)

    else
        Sigma <- matrix(c(param["cov11"], param["cov12"], param["cov13"],
                          param["cov12"], param["cov22"], param["cov23"],
                          param["cov13"], param["cov23"], param["cov33"]),
                        3, 3)

    iSigma <- solve(Sigma)

    if (iso)
        ext.coeff <- function(h)
            2 * pnorm(sqrt(h^2 / param["cov11"]) / 2)

    else
        ext.coeff <- function(h)
            2 * pnorm(sqrt(h %*% iSigma %*% h) / 2)

    fitted <- list(fitted.values = opt$par, std.err = std.err,
                   var.cov = var.cov, fixed = unlist(fixed.param), param = param,
                   deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                   counts = opt$counts, message = opt$message, data = data, est = "MPLE",
                   logLik = -opt$value, opt.value = opt$value, model = "Smith",
                   fit.marge = fit.marge, ext.coeff = ext.coeff, cov.mod = "Gaussian",
                   lik.fun = nllh, coord = coord, ihessian = ihessian, var.score = var.score,
                   marg.cov = NULL, nllh = nllh, iso = iso, weighted = weighted)

    class(fitted) <- c(fitted$model, "maxstab")
    return(fitted)
}

##This functions fits the model from a generic R formula
##i.e. classical linear models as well as smoothing splines with
##radial basis functions may be used to model spatially the GEV
##parameters
smithform <- function(data, coord, loc.form, scale.form, shape.form, start, fit.marge = TRUE,
                      iso = TRUE, marg.cov = NULL, ..., warn = TRUE, method = "BFGS",
                      control = list(), corr = FALSE,
                      weights = weights, temp.cov = NULL, temp.form.loc = NULL,
                      temp.form.scale = NULL, temp.form.shape = NULL, check.grad = FALSE){
    ##data is a matrix with each column corresponds to one location
    ##coord is a matrix giving the coordinates (1 row = 1 station)
    n.site <- ncol(data)
    n.obs <- nrow(data)
    dist.dim <- ncol(coord)
    n.pair <- n.site * (n.site - 1) / 2

    distVec <- distance(coord, vec = TRUE)
    weighted <- !is.null(weights)

    if (!weighted)
        ##Set the weights to 0 as they won't be used anyway
        weights <- 0

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

    if (is.null(marg.cov))
        covariables <- data.frame(coord)

    else
        covariables <- data.frame(coord, marg.cov)

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

    loc.type <- loc.model$type
    scale.type <- scale.model$type
    shape.type <- shape.model$type

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

    if (dist.dim == 2){
        if (iso)
            param <- c("cov", param)

        else
            param <- c("cov11", "cov12", "cov22", param)
    }

    else{
        if (iso)
            param <- c("cov", param)

        else
            param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33", param)
    }

    ##First create a "void" function
    nplk <- function(x) x

    ##And define the "body" of the function as the number of parameters
    ##to estimate depends on n.site
    if (dist.dim == 2){

        if (iso)
            body(nplk) <- parse(text = paste("-.C('smithdsgnmat', as.double(data), as.double(distVec),
as.integer(n.site), as.integer(n.obs), as.integer(weighted), as.double(weights),
as.double(loc.dsgn.mat), as.double(loc.pen.mat), as.integer(n.loccoeff), as.integer(n.pparloc),
as.double(loc.penalty), as.double(scale.dsgn.mat), as.double(scale.pen.mat),
as.integer(n.scalecoeff), as.integer(n.pparscale), as.double(scale.penalty),
as.double(shape.dsgn.mat), as.double(shape.pen.mat), as.integer(n.shapecoeff),
as.integer(n.pparshape), as.double(shape.penalty), as.integer(use.temp.cov),
as.double(temp.dsgn.mat.loc), as.double(temp.pen.mat.loc), as.integer(n.tempcoeff.loc),
as.integer(n.ppartemp.loc), as.double(temp.penalty.loc), as.double(temp.dsgn.mat.scale),
as.double(temp.pen.mat.scale), as.integer(n.tempcoeff.scale), as.integer(n.ppartemp.scale),
as.double(temp.penalty.scale), as.double(temp.dsgn.mat.shape), as.double(temp.pen.mat.shape),
as.integer(n.tempcoeff.shape), as.integer(n.ppartemp.shape), as.double(temp.penalty.shape),",
                                paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.loc, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.scale, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.shape, collapse = ","), ")), "),
                                "as.double(cov), as.double(0.0), as.double(cov), dns = double(1),
PACKAGE = 'SpatialExtremes', NAOK = TRUE)$dns"))

        else
            body(nplk) <- parse(text = paste("-.C('smithdsgnmat', as.double(data), as.double(distVec),
as.integer(n.site), as.integer(n.obs), as.integer(weighted), as.double(weights),
as.double(loc.dsgn.mat), as.double(loc.pen.mat), as.integer(n.loccoeff), as.integer(n.pparloc),
as.double(loc.penalty), as.double(scale.dsgn.mat), as.double(scale.pen.mat),
as.integer(n.scalecoeff), as.integer(n.pparscale), as.double(scale.penalty),
as.double(shape.dsgn.mat), as.double(shape.pen.mat), as.integer(n.shapecoeff),
as.integer(n.pparshape), as.double(shape.penalty), as.integer(use.temp.cov),
as.double(temp.dsgn.mat.loc), as.double(temp.pen.mat.loc), as.integer(n.tempcoeff.loc),
as.integer(n.ppartemp.loc), as.double(temp.penalty.loc), as.double(temp.dsgn.mat.scale),
as.double(temp.pen.mat.scale), as.integer(n.tempcoeff.scale), as.integer(n.ppartemp.scale),
as.double(temp.penalty.scale), as.double(temp.dsgn.mat.shape), as.double(temp.pen.mat.shape),
as.integer(n.tempcoeff.shape), as.integer(n.ppartemp.shape), as.double(temp.penalty.shape),",
                                paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.loc, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.scale, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.shape, collapse = ","), ")), "),
                                "as.double(cov11), as.double(cov12), as.double(cov22), dns = double(1),
PACKAGE = 'SpatialExtremes', NAOK = TRUE)$dns"))
    }

    else{
        if (iso)
            body(nplk) <- parse(text = paste("-.C('smithdsgnmat3d', as.double(data), as.double(distVec),
as.integer(n.site), as.integer(n.obs), as.integer(weighted), as.double(weights),
as.double(loc.dsgn.mat), as.double(loc.pen.mat), as.integer(n.loccoeff), as.integer(n.pparloc),
as.double(loc.penalty), as.double(scale.dsgn.mat), as.double(scale.pen.mat),
as.integer(n.scalecoeff), as.integer(n.pparscale), as.double(scale.penalty),
as.double(shape.dsgn.mat), as.double(shape.pen.mat), as.integer(n.shapecoeff),
as.integer(n.pparshape), as.double(shape.penalty), as.integer(use.temp.cov),
as.double(temp.dsgn.mat.loc), as.double(temp.pen.mat.loc), as.integer(n.tempcoeff.loc),
as.integer(n.ppartemp.loc), as.double(temp.penalty.loc), as.double(temp.dsgn.mat.scale),
as.double(temp.pen.mat.scale), as.integer(n.tempcoeff.scale), as.integer(n.ppartemp.scale),
as.double(temp.penalty.scale), as.double(temp.dsgn.mat.shape), as.double(temp.pen.mat.shape),
as.integer(n.tempcoeff.shape), as.integer(n.ppartemp.shape), as.double(temp.penalty.shape),",
                                paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.loc, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.scale, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.shape, collapse = ","), ")), "),
                                "as.double(cov), as.double(0.0), as.double(0.0), as.double(cov),
as.double(0.0), as.double(cov), dns = double(1), PACKAGE = 'SpatialExtremes', NAOK = TRUE)$dns"))

        else
            body(nplk) <- parse(text = paste("-.C('smithdsgnmat3d', as.double(data), as.double(distVec),
as.integer(n.site), as.integer(n.obs), as.integer(weighted), as.double(weights),
as.double(loc.dsgn.mat), as.double(loc.pen.mat), as.integer(n.loccoeff), as.integer(n.pparloc),
as.double(loc.penalty), as.double(scale.dsgn.mat), as.double(scale.pen.mat),
as.integer(n.scalecoeff), as.integer(n.pparscale), as.double(scale.penalty),
as.double(shape.dsgn.mat), as.double(shape.pen.mat), as.integer(n.shapecoeff),
as.integer(n.pparshape), as.double(shape.penalty), as.integer(use.temp.cov),
as.double(temp.dsgn.mat.loc), as.double(temp.pen.mat.loc), as.integer(n.tempcoeff.loc),
as.integer(n.ppartemp.loc), as.double(temp.penalty.loc), as.double(temp.dsgn.mat.scale),
as.double(temp.pen.mat.scale), as.integer(n.tempcoeff.scale), as.integer(n.ppartemp.scale),
as.double(temp.penalty.scale), as.double(temp.dsgn.mat.shape), as.double(temp.pen.mat.shape),
as.integer(n.tempcoeff.shape), as.integer(n.ppartemp.shape), as.double(temp.penalty.shape),",
                                paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.loc, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.scale, collapse = ","), ")), "),
                                paste("as.double(c(", paste(temp.names.shape, collapse = ","), ")), "),
                                "as.double(cov11), as.double(cov12), as.double(cov13),
as.double(cov22), as.double(cov23), as.double(cov33), dns = double(1),
PACKAGE = 'SpatialExtremes', NAOK = TRUE)$dns"))
    }

    ##Define the formal arguments of the function
    form.nplk <- NULL
    for (i in 1:length(param))
        form.nplk <- c(form.nplk, alist(a=))

    names(form.nplk) <- param
    formals(nplk) <- form.nplk

    if (missing(start)) {

        start <- .start.smith(data, coord, covariables, loc.form, scale.form,
                              shape.form, method = method, iso = iso,
                              ...)

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

        start <- c(start, as.list(c(tempCoeff.loc, tempCoeff.scale, tempCoeff.shape)))
        start <- start[!(param %in% names(list(...)))]

    }

    if (!is.list(start))
        stop("'start' must be a named list")

    if (!length(start))
        stop("there are no parameters left to maximize over")

    nm <- names(start)
    l <- length(nm)
    f <- formals(nplk)
    names(f) <- param
    m <- match(nm, param)

    if(any(is.na(m)))
        stop("'start' specifies unknown arguments")

    formals(nplk) <- c(f[m], f[-m])
    nllh <- function(p, ...) nplk(p, ...)

    if(l > 1)
        body(nllh) <- parse(text = paste("nplk(", paste("p[", 1:l, "]", collapse = ", "), ", ...)"))

    fixed.param <- list(...)[names(list(...)) %in% param]

    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")

    start.arg <- c(list(p = unlist(start)), fixed.param)

    init.lik <- do.call("nllh", start.arg)
    if (warn && (init.lik >= 1.0e15))
        warning("negative log-likelihood is infinite at starting values")

    if (method == "nlminb"){
        start <- as.numeric(start)
        opt <- nlminb(start, nllh, ..., control = control)
        opt$counts <- opt$evaluations
        opt$value <- opt$objective
        names(opt$par) <- nm

        if ((opt$convergence != 0) || (opt$value >= 1.0e15)) {
            if (warn)
                warning("optimization may not have succeeded")
        }

        if (opt$convergence == 0)
            opt$convergence <- "successful"
    }

    if (method == "nlm"){
        start <- as.numeric(start)
        opt <- nlm(nllh, start, ...)
        opt$counts <- opt$iterations
        names(opt$counts) <- "function"
        opt$value <- opt$minimum
        opt$par <- opt$estimate
        names(opt$par) <- nm

        if (opt$code <= 2)
            opt$convergence <- "sucessful"

        if (opt$code == 3)
            opt$convergence <- "local minimum or 'steptol' is too small"

        if (opt$code == 4)
            opt$convergence <- "iteration limit reached"

        if (opt$code == 5)
            opt$convergence <- "optimization failed"
    }

    if(!(method %in% c("nlm", "nlminb"))) {

        opt <- optim(start, nllh, ..., method = method, control = control)

        if ((opt$convergence != 0) || (opt$value >= 1.0e15)) {
            if (warn)
                warning("optimization may not have succeeded")

            if (opt$convergence != 0)
                opt$convergence <- "iteration limit reached"
        }

        else opt$convergence <- "successful"
    }

    if (opt$value == init.lik){
        if (warn)
            warning("optimization stayed at the starting values.")

        opt$convergence <- "Stayed at start. val."
    }

    param.names <- param
    param <- c(opt$par, unlist(fixed.param))
    param <- param[param.names]

    if (iso){
        if (dist.dim == 2){
            param <- c(param["cov"], 0, param["cov"], param[-1])
            names(param)[1:3] <- c("cov11", "cov12", "cov22")
        }

        else{
            param <- c(param["cov"], 0, 0, param["cov"], 0, param["cov"],
                       param[-1])
            names(param)[1:6] <- c("cov11", "cov12", "cov13", "cov22", "cov23",
                                   "cov33")
        }
    }

    ##Reset the weights to their original values
    if ((length(weights) == 1) && (weights == 0))
        weights <- NULL

    std.err <- .smithstderr(param, data, distVec, loc.dsgn.mat, scale.dsgn.mat, shape.dsgn.mat,
                            temp.dsgn.mat.loc, temp.dsgn.mat.scale, temp.dsgn.mat.shape, use.temp.cov,
                            fit.marge = fit.marge, fixed.param = names(fixed.param),
                            param.names = param.names, iso = iso, weights = weights)

    if (check.grad)
        print(round(rbind(numerical = -opt$grad, analytical = std.err$grad), 3))

    opt$hessian <- std.err$hessian
    var.score <- std.err$var.score
    ihessian <- try(solve(opt$hessian), silent = TRUE)

    if(!is.matrix(ihessian) || any(is.na(var.score))){
        if (warn)
            warning("Observed information matrix is singular. No standard error will be computed.")

        std.err.type <- "none"
        return
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

        colnames(var.cov) <- rownames(var.cov) <- colnames(ihessian) <-
            rownames(ihessian) <- names(std.err) <- nm
    }

    if (std.err.type == "none"){
        std.err <- std.err.type <- corr.mat <- NULL
        var.cov <- ihessian <- var.score <- NULL
    }

    if (dist.dim == 2)
        Sigma <- matrix(c(param["cov11"], param["cov12"], param["cov12"],
                          param["cov22"]), 2, 2)

    else
        Sigma <- matrix(c(param["cov11"], param["cov12"], param["cov13"],
                          param["cov12"], param["cov22"], param["cov23"],
                          param["cov13"], param["cov23"], param["cov33"]),
                        3, 3)

    iSigma <- solve(Sigma)

    if (iso)
        ext.coeff <- function(h)
            2 * pnorm(sqrt(h^2 / param["cov11"]) / 2)
    else
        ext.coeff <- function(h)
            2 * pnorm(sqrt(h %*% iSigma %*% h) / 2)

    fitted <- list(fitted.values = opt$par, std.err = std.err,
                   var.cov = var.cov, fixed = unlist(fixed.param), param = param,
                   deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                   counts = opt$counts, message = opt$message, data = data, est = "MPLE",
                   logLik = -opt$value, opt.value = opt$value, model = "Smith", coord = coord,
                   fit.marge = fit.marge, ext.coeff = ext.coeff, cov.mod = "Gaussian",
                   loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                   lik.fun = nllh, loc.type = loc.type, scale.type = scale.type,
                   shape.type = shape.type, ihessian = ihessian, var.score = var.score,
                   marg.cov = marg.cov, nllh = nllh, iso = iso, weighted = weighted)

    class(fitted) <- c(fitted$model, "maxstab")
    return(fitted)
}
