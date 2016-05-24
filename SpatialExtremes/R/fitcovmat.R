lsmaxstab <- function(data, coord, cov.mod = "gauss", marge = "emp", control = list(),
                      iso = FALSE, ..., weighted = TRUE){

    if (cov.mod == "gauss")
        fitcovmat(data, coord, marge, control = control, ..., iso = iso,
                  weighted = weighted)

    else
        fitcovariance(data, coord, cov.mod, marge = marge, control = control,
                      ..., weighted = weighted)
}

fitcovmat <- function(data, coord, marge = "emp", iso = FALSE, control = list(),
                      ..., start, weighted = TRUE){

    if (is.null(control$maxit))
        control$maxit <- 10000

    n.site <- ncol(data)
    n.pairs <- n.site * (n.site - 1) / 2
    dist.dim <- ncol(coord)

    extcoeff <- fitextcoeff(data, coord, estim = "Smith", plot = FALSE, loess = FALSE,
                            marge = marge, std.err = weighted)

    if (weighted){
        weights <- extcoeff[,"std.err"]
        ##Check if there are really small weights, this could happen in few
        ##cases
        weights[weights <= 1e-4] <- mean(weights)
    }

    else
        weights <- rep(1, n.pairs)

    extcoeff <- extcoeff[,"ext.coeff"]
    dist <- distance(coord, vec = TRUE)

    if (dist.dim == 2){

        if (iso){
            param <- "cov"

            fun2diso <- function(cov)
                .C("fitcovmat2d", as.double(cov), as.double(0.0),
                   as.double(cov), as.integer(n.pairs), as.double(dist),
                   as.double(extcoeff), as.double(weights), ans = double(1),
                   PACKAGE = "SpatialExtremes")$ans
        }

        else{
            param <- c("cov11", "cov12", "cov22")

            fun2d <- function(cov11, cov12, cov22)
                .C("fitcovmat2d", as.double(cov11), as.double(cov12),
                   as.double(cov22), as.integer(n.pairs), as.double(dist),
                   as.double(extcoeff), as.double(weights), ans = double(1),
                   PACKAGE = "SpatialExtremes")$ans
        }
    }

    if (dist.dim == 3){

        if (iso){
            param <- "cov"

            fun3diso <- function(cov)
                .C("fitcovmat3d", as.double(cov), as.double(0.0), as.double(0.0),
                   as.double(cov), as.double(0.0), as.double(cov),
                   as.integer(n.pairs), as.double(dist), as.double(extcoeff),
                   as.double(weights), ans = double(1), PACKAGE = "SpatialExtremes")$ans
        }

        else{
            param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33")

            fun3d <- function(cov11, cov12, cov13, cov22, cov23, cov33)
                .C("fitcovmat3d", as.double(cov11), as.double(cov12), as.double(cov13),
                   as.double(cov22), as.double(cov23), as.double(cov33),
                   as.integer(n.pairs), as.double(dist), as.double(extcoeff),
                   as.double(weights), ans = double(1), PACKAGE = "SpatialExtremes")$ans
        }
    }

    fixed.param <- list(...)[names(list(...)) %in% param]

    if (missing(start)){

        a <- 4 * qnorm(pmin(extcoeff, 2) / 2)^2
        sigma.start <- mean(rowSums(dist^2) / a)

        if (iso)
            start <- list(cov = sigma.start)

        else{
            if (dist.dim == 2){
                if (any(names(fixed.param) == "cov12"))
                    start <- list(cov11 = 1 + 2 * abs(list(...)$cov12),
                                  cov12 = list(...)$cov12,
                                  cov22 = 1 + 2 * abs(list(...)$cov12))

                else
                    start <- list(cov11 = sigma.start, cov12 = 0, cov22 = sigma.start)
            }

            if (dist.dim == 3)
                start <- list(cov11 = sigma.start, cov12 = 0, cov13 = 0, cov22 = sigma.start,
                              cov23 = 0, cov33 = sigma.start)
        }

        start <- start[!(param %in% names(list(...)))]
    }

    if (!is.list(start))
        stop("'start' must be a named list")

    if (!length(start))
        stop("there are no parameters left to maximize over")

    nm <- names(start)
    l <- length(nm)

    if (dist.dim == 2){
        if (iso)
            f <- formals(fun2diso)

        else
            f <- formals(fun2d)
    }

    if (dist.dim == 3){
        if (iso)
            f <- formals(fun3diso)

        else
            f <- formals(fun3d)
    }

    names(f) <- param
    m <- match(nm, param)

    if(any(is.na(m)))
        stop("'start' specifies unknown arguments")

    if (dist.dim == 2){
        if (iso){
            formals(fun2diso) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) fun2diso(p, ...)


            if (l > 1)
                body(obj.fun) <- parse(text = paste("fun2diso(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }

        else{
            formals(fun2d) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) fun2d(p, ...)


            if (l > 1)
                body(obj.fun) <- parse(text = paste("fun2d(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }
    }

    if (dist.dim == 3){
        if (iso){
            formals(fun3diso) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) fun3diso(p, ...)


            if (l > 1)
                body(obj.fun) <- parse(text = paste("fun3diso(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }

        else{
            formals(fun3d) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) fun3d(p, ...)


            if (l > 1)
                body(obj.fun) <- parse(text = paste("fun3d(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }
    }

    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")

    opt <- optim(unlist(start), obj.fun, hessian = FALSE, control = control,
                 ...)

    if (opt$convergence == 0)
        opt$convergence <- "successful"

    else if (opt$convergence == 1)
        opt$convergence <- "iteration limit reached"

    else
        opt$convergence <- "Optimization may have failed"

    param.names <- param
    names(opt$par) <- names(start)
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
            2 * pnorm(0.5 * h / sqrt(param["cov11"]))

    else
        ext.coeff <- function(h)
            2 * pnorm(0.5 * sqrt(h %*% iSigma %*% h))

    fitted <- list(fitted.values = opt$par, fixed = unlist(fixed.param),
                   param = param, convergence = opt$convergence,
                   counts = opt$counts, message = opt$message, data = data,
                   est = "Least Squares", opt.value = opt$value, model = "Smith",
                   coord = coord, fit.marge = FALSE, cov.mod = "Gaussian",
                   ext.coeff = ext.coeff, iso = iso, weighted = weighted)

    class(fitted) <- c(fitted$model, "maxstab")
    return(fitted)

}

fitcovariance <- function(data, coord, cov.mod, marge = "emp", control = list(),
                          ..., start, weighted = TRUE){

    if (is.null(control$maxit))
        control$maxit <- 10000

    n.site <- ncol(data)
    n.pairs <- n.site * (n.site - 1) / 2
    dist.dim <- ncol(coord)

    if (substr(cov.mod, 1, 1) %in% c("i", "g", "t")){

        if (substr(cov.mod, 1, 1) == "i")
            model <- "iSchlather"

        else if (substr(cov.mod, 1, 1) == "g")
            model <- "Geometric"

        else
            model <- "Extremal-t"

        cov.mod <- substr(cov.mod, 2, 8)
    }

    else if (cov.mod == "brown")
        model <- "Brown-Resnick"

    else
        model <- "Schlather"

    if (!(cov.mod %in% c("whitmat","cauchy","powexp","bessel","caugen", "brown")))
        stop("''cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp', 'bessel', 'caugen', 'brown'")

    if (cov.mod == "whitmat")
        cov.mod.num <- 1
    if (cov.mod == "cauchy")
        cov.mod.num <- 2
    if (cov.mod == "powexp")
        cov.mod.num <- 3
    if (cov.mod == "bessel")
        cov.mod.num <- 4
    if (cov.mod == "caugen")
        cov.mod.num <- 5

    extcoeff <- fitextcoeff(data, coord, estim = "Smith", plot = FALSE, loess = FALSE,
                            marge = marge, std.err = weighted)

    if (weighted){
        weights <- extcoeff[,"std.err"]
        ##Check if there are really small weights, this could happen in few
        ##cases
        weights[weights <= 1e-4] <- mean(weights)
    }

    else
        weights <- rep(1, n.pairs)

    extcoeff <- extcoeff[,"ext.coeff"]
    dist <- distance(coord)

    if (cov.mod == "brown")
        param <- c("range", "smooth")

    else
        param <- c("nugget", "range", "smooth")

    if (model == "iSchlather")
        param <- c("alpha", param)

    if (model == "Geometric"){
        param <- c("sigma2", param)

        if (is.null(control$sigma2Bound))
            sigma2Bound <- 50

        else
            sigma2Bound <- control$sigma2Bound
    }

    if (cov.mod == "caugen")
        param <- c(param, "smooth2")

    if (model == "Extremal-t")
        param <- c(param, "DoF")

    if (model == "Schlather"){
        if (cov.mod == "caugen")
            funS2 <- function(nugget, range, smooth, smooth2)
                .C("fitcovariance", as.integer(cov.mod.num), as.double(nugget), as.double(range),
                   as.double(smooth), as.double(smooth2), as.integer(n.pairs), as.integer(dist.dim),
                   as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1),
                   PACKAGE = "SpatialExtremes")$ans

        else
            funS1 <- function(nugget, range, smooth)
                .C("fitcovariance", as.integer(cov.mod.num), as.double(nugget), as.double(range),
                   as.double(smooth), double(1), as.integer(n.pairs), as.integer(dist.dim),
                   as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1),
                   PACKAGE = "SpatialExtremes")$ans
    }

    else if (model == "iSchlather"){
        if (cov.mod == "caugen")
            funI2 <- function(alpha, nugget, range, smooth, smooth2)
                .C("fiticovariance", as.integer(cov.mod.num), as.double(alpha), as.double(nugget),
                   as.double(range), as.double(smooth), as.double(smooth2), as.integer(n.pairs),
                   as.integer(dist.dim), as.double(dist), as.double(extcoeff), as.double(weights),
                   ans = double(1), PACKAGE = "SpatialExtremes")$ans

        else
            funI1 <- function(alpha, nugget, range, smooth)
                .C("fiticovariance", as.integer(cov.mod.num), as.double(alpha), as.double(nugget),
                   as.double(range), as.double(smooth), double(1), as.integer(n.pairs), as.integer(dist.dim),
                   as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1),
                   PACKAGE = "SpatialExtremes")$ans
    }

    else if (model == "Geometric"){
        if (cov.mod == "caugen")
            funG2 <- function(sigma2, nugget, range, smooth, smooth2)
                .C("fitgcovariance", as.integer(cov.mod.num), as.double(sigma2), as.double(sigma2Bound),
                   as.double(nugget), as.double(range), as.double(smooth), as.double(smooth2),
                   as.integer(n.pairs), as.integer(dist.dim), as.double(dist), as.double(extcoeff),
                   as.double(weights), ans = double(1), PACKAGE = "SpatialExtremes")$ans

        else
            funG1 <- function(sigma2, nugget, range, smooth)
                .C("fitgcovariance", as.integer(cov.mod.num), as.double(sigma2), as.double(sigma2Bound),
                   as.double(nugget), as.double(range), as.double(smooth), double(1), as.integer(n.pairs),
                   as.integer(dist.dim), as.double(dist), as.double(extcoeff), as.double(weights),
                   ans = double(1), PACKAGE = "SpatialExtremes")$ans
    }

    else if (model == "Extremal-t") {
        if (cov.mod == "caugen")
            funT2 <- function(nugget, range, smooth, smooth2, DoF)
                .C("fittcovariance", as.integer(cov.mod.num), as.double(nugget), as.double(range), as.double(smooth),
                   as.double(smooth2), as.double(DoF), as.integer(n.pairs), as.integer(dist.dim), as.double(dist),
                   as.double(extcoeff), as.double(weights), ans = double(1), PACKAGE = "SpatialExtremes")$ans

        else
            funT1 <- function(nugget, range, smooth, DoF)
                .C("fittcovariance", as.integer(cov.mod.num), as.double(nugget), as.double(range), as.double(smooth),
                   double(1), as.double(DoF), as.integer(n.pairs), as.integer(dist.dim), as.double(dist),
                   as.double(extcoeff), as.double(weights), ans = double(1), PACKAGE = "SpatialExtremes")$ans
    }

    else
        funBR <- function(range, smooth)
            .C("fitbrcovariance", as.double(range), as.double(smooth), as.integer(n.pairs),
               as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1),
               PACKAGE = "SpatialExtremes")$ans

    fixed.param <- list(...)[names(list(...)) %in% param]

    if (missing(start)){

        if (cov.mod == "whitmat"){
            range.start <- max(dist) / 2.33
            smooth.start <- 0.5
        }

        if (cov.mod == "cauchy"){
            range.start <- max(dist) / 10
            smooth.start <- 0.35
        }

        if (cov.mod == "bessel"){
            range.start <- 0.75 * max(dist) / 8
            smooth.start <- 2
        }

        if (cov.mod == "caugen"){
            range.start <- 0.75 * max(dist) / 20
            smooth.start <- 1
        }

        if (cov.mod == "powexp"){
            range.start <- max(dist) / 2.33
            smooth.start <- 1
        }

        if (cov.mod == "brown"){
            range.start <- max(dist) / 2.73
            smooth.start <- 0.75
        }

        if (cov.mod == "brown")
            start <- list(range = range.start, smooth = smooth.start)

        else
            start <- list(nugget = .1, range = range.start, smooth = smooth.start)

        if (cov.mod == "caugen")
            start <- c(start, list(smooth2 = 0.5))

        if (model == "iSchlather")
            start <- c(list(alpha = 0.5), start)

        if (model == "Geometric")
            start <- c(list(sigma2 = 9), start)

        if (model == "Extremal-t"){
          start <- list(nugget = 0.1, range = max(dist) / 2.33, smooth = 1, DoF = 1)

          if (cov.mod == "caugen")
            start <- list(nugget = 0.1, range = 2 * max(dist), smooth = 1, smooth2 = 1, DoF = 1)

          if (cov.mod == "cauchy")
            start <- list(nugget = 0.1, range = 0.07 * max(dist), smooth = 0.1, DoF = 1)
        }

        start <- start[!(param %in% names(list(...)))]
    }

    if (!is.list(start))
        stop("'start' must be a named list")

    if (!length(start))
        stop("there are no parameters left to maximize over")

    nm <- names(start)
    l <- length(nm)

    if (model == "Schlather"){
        if (cov.mod == "caugen")
            f <- formals(funS2)

        else
            f <- formals(funS1)
    }

    else if (model == "iSchlather"){
        if (cov.mod == "caugen")
            f <- formals(funI2)

        else
            f <- formals(funI1)
    }

    else if (model == "Geometric"){
        if (cov.mod == "caugen")
            f <- formals(funG2)

        else
            f <- formals(funG1)
    }

    else if (model == "Extremal-t"){
        if (cov.mod == "caugen")
            f <- formals(funT2)

        else
            f <- formals(funT1)
    }

    else
        f <- formals(funBR)

    names(f) <- param
    m <- match(nm, param)

    if(any(is.na(m)))
        stop("'start' specifies unknown arguments")

    if (model == "Schlather"){
        if (cov.mod == "caugen"){
            formals(funS2) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) funS2(p, ...)

            if (l > 1)
                body(obj.fun) <- parse(text = paste("funS2(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }

        else {
            formals(funS1) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) funS1(p, ...)

            if (l > 1)
                body(obj.fun) <- parse(text = paste("funS1(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }
    }

    else if (model == "iSchlather"){
        if (cov.mod == "caugen"){
            formals(funI2) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) funI2(p, ...)

            if (l > 1)
                body(obj.fun) <- parse(text = paste("funI2(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }

        else {
            formals(funI1) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) funI1(p, ...)

            if (l > 1)
                body(obj.fun) <- parse(text = paste("funI1(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }
    }

    else if (model == "Geometric"){

        if (cov.mod == "caugen"){
            formals(funG2) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) funG2(p, ...)

            if (l > 1)
                body(obj.fun) <- parse(text = paste("funG2(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }

        else {
            formals(funG1) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) funG1(p, ...)

            if (l > 1)
                body(obj.fun) <- parse(text = paste("funG1(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }
    }

    else if (model == "Extremal-t"){
        if (cov.mod == "caugen"){
            formals(funT2) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) funT2(p, ...)

            if (l > 1)
                body(obj.fun) <- parse(text = paste("funT2(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }

        else {
            formals(funT1) <- c(f[m], f[-m])
            obj.fun <- function(p, ...) funT1(p, ...)

            if (l > 1)
                body(obj.fun) <- parse(text = paste("funT1(", paste("p[",1:l,
                                       "]", collapse = ", "), ", ...)"))
        }
    }


    else {
        formals(funBR) <- c(f[m], f[-m])
        obj.fun <- function(p, ...) funBR(p, ...)

        if (l > 1)
            body(obj.fun) <- parse(text = paste("funBR(", paste("p[",1:l,
                                   "]", collapse = ", "), ", ...)"))
    }

    if(any(!(param %in% c(nm,names(fixed.param)))))
        stop("unspecified parameters")

    opt <- optim(unlist(start), obj.fun, hessian = FALSE, control = control, ...)

    if (opt$convergence == 1)
        opt$convergence <- "iteration limit reached"

    else if (opt$convergence == 0)
        opt$convergence <- "successful"

    else
        opt$convergence <- "Optimization may have failed"

    param.names <- param
    names(opt$par) <- names(start)
    param <- c(opt$par, unlist(fixed.param))
    param <- param[param.names]

    if (cov.mod != "brown"){
        if (cov.mod != "caugen")
            cov.fun <- covariance(nugget = param["nugget"], sill = 1 - param["nugget"], range = param["range"],
                                  smooth = param["smooth"], cov.mod = cov.mod, plot = FALSE)

        else
            cov.fun <- covariance(nugget = param["nugget"], sill = 1 - param["nugget"], range = param["range"],
                                  smooth = param["smooth"], smooth2 = param["smooth2"],
                                  cov.mod = cov.mod, plot = FALSE)

        if (model == "Schlather")
            ext.coeff <- function(h)
                1 + sqrt(0.5 - 0.5 * cov.fun(h))

        else if (model == "iSchlather")
            ext.coeff <- function(h)
                2 * param["alpha"] + (1 - param["alpha"]) * (1 + sqrt(0.5 - 0.5 * cov.fun(h)))

        else if (model == "Extremal-t")
            ext.coeff <- function(h)
                2 * pt(sqrt((1 - cov.fun(h)) * (param["DoF"] + 1) / (1 + cov.fun(h))), param["DoF"] + 1)

        else
            ext.coeff <- function(h)
                2 * pnorm(sqrt(0.5 * param["sigma2"] * (1 - cov.fun(h))))
    }

    else {
        cov.fun <- NA
        ext.coeff <- function(h)
            2 * pnorm((h / param["range"])^(0.5 * param["smooth"]) / sqrt(2))
    }

    fitted <- list(fitted.values = opt$par, fixed = unlist(fixed.param),
                   param = param, convergence = opt$convergence,
                   counts = opt$counts, message = opt$message, data = data,
                   est = "Least Squares", opt.value = opt$value, model = model,
                   coord = coord, fit.marge = FALSE, cov.mod = cov.mod,
                   cov.fun = cov.fun, ext.coeff = ext.coeff, weighted = weighted)

    class(fitted) <- c(fitted$model, "maxstab")
    return(fitted)

}
