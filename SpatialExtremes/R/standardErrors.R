.smithstderr <- function(par, data, distVec, loc.dsgn.mat, scale.dsgn.mat, shape.dsgn.mat,
                         temp.dsgn.mat.loc, temp.dsgn.mat.scale, temp.dsgn.mat.shape,
                         use.temp.cov, fit.marge, fixed.param, param.names,
                         iso = TRUE, weights){

    ##data is a matrix with each column corresponds to one location
    ##distVec is the a matrix giving the "distance vector" for each pair
    ##(1 row = 1 station)
    n.site <- ncol(data)
    n.obs <- nrow(data)
    n.pairs <- n.site * (n.site - 1) / 2
    dist.dim <- ncol(distVec)
    n.param <- length(param.names)

    if (is.null(weights))
        weights <- rep(1, n.pairs)

    if (iso){
        if (dist.dim == 2)
            n.param <- n.param + 2

        else
            n.param <- n.param + 5
    }

    cov11 <- par["cov11"]
    cov12 <- par["cov12"]
    cov22 <- par["cov22"]

    if (dist.dim == 3){
        cov13 <- par["cov13"]
        cov23 <- par["cov23"]
        cov33 <- par["cov33"]
    }

    if (fit.marge){

        n.loccoeff <- ncol(loc.dsgn.mat)
        n.scalecoeff <- ncol(scale.dsgn.mat)
        n.shapecoeff <- ncol(shape.dsgn.mat)

        loc.idx <- which(substr(names(par), 1, 3) == "loc")
        scale.idx <- which(substr(names(par), 1, 5) == "scale")
        shape.idx <- which(substr(names(par), 1, 5) == "shape")

        loc.param <- par[loc.idx]
        scale.param <- par[scale.idx]
        shape.param <- par[shape.idx]

        n.temploccoeff <- ifelse(use.temp.cov[1], ncol(temp.dsgn.mat.loc), 0)
        n.tempscalecoeff <- ifelse(use.temp.cov[2], ncol(temp.dsgn.mat.scale), 0)
        n.tempshapecoeff <- ifelse(use.temp.cov[3], ncol(temp.dsgn.mat.scale), 0)

        temp.loc.idx <- which(substr(names(par), 1, 12) == "tempCoeffLoc")
        temp.scale.idx <- which(substr(names(par), 1, 14) == "tempCoeffScale")
        temp.shape.idx <- which(substr(names(par), 1, 14) == "tempCoeffShape")

        temp.loc.param <- ifelse(length(temp.loc.idx) > 0, par[temp.loc.idx], 0)
        temp.scale.param <- ifelse(length(temp.scale.idx) > 0, par[temp.scale.idx], 0)
        temp.shape.param <- ifelse(length(temp.shape.idx) > 0, par[temp.shape.idx], 0)
    }

    else {
        n.loccoeff <- n.scalecoeff <- n.shapecoeff <- loc.param <- scale.param <- shape.param <- 1
        temp.loc.param <- temp.scale.param <- temp.shape.param <- n.temploccoeff <-
            n.tempscalecoeff <- n.tempshapecoeff <- 0
    }

    if (dist.dim == 2)
        std.err <- .C("smithstderr", as.double(data), as.double(distVec), as.integer(n.site),
                      as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
                      as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                      as.integer(n.shapecoeff), as.double(temp.dsgn.mat.loc), as.integer(n.temploccoeff),
                      as.double(temp.dsgn.mat.scale), as.integer(n.tempscalecoeff),
                      as.double(temp.dsgn.mat.shape), as.integer(n.tempshapecoeff),
                      as.double(loc.param), as.double(scale.param), as.double(shape.param),
                      as.double(temp.loc.param), as.double(temp.scale.param), as.double(temp.shape.param),
                      as.double(cov11), as.double(cov12), as.double(cov22), as.integer(fit.marge),
                      as.integer(use.temp.cov), as.double(weights), hess = double(n.obs * n.param * n.pairs),
                      grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes", NAOK = TRUE)

    else
        std.err <- .C("smithstderr3d", as.double(data), as.double(distVec), as.integer(n.site),
                      as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
                      as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                      as.integer(n.shapecoeff), as.double(temp.dsgn.mat.loc), as.integer(n.temploccoeff),
                      as.double(temp.dsgn.mat.scale), as.integer(n.tempscalecoeff), as.double(temp.dsgn.mat.shape),
                      as.integer(n.tempshapecoeff), as.double(loc.param), as.double(scale.param),
                      as.double(shape.param), as.double(temp.loc.param), as.double(temp.scale.param),
                      as.double(temp.shape.param), as.double(cov11), as.double(cov12), as.double(cov13),
                      as.double(cov22), as.double(cov23), as.double(cov33), fit.marge, as.integer(use.temp.cov),
                      as.double(weights), hess = double(n.obs * n.param * n.pairs), grad = double(n.obs * n.param),
                      PACKAGE = "SpatialExtremes", NAOK = TRUE)

    grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
    hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)

    if (iso){
        if (dist.dim == 2){
            grad[,1] <- rowSums(grad[,c(1,3)])
            grad <- grad[,-(2:3), drop = FALSE]

            hess[,1] <- rowSums(hess[,c(1,3)])
            hess <- hess[,-(2:3), drop = FALSE]
        }

        if (dist.dim == 3){
            grad[,1] <- rowSums(grad[,c(1,4,6)])
            grad <- grad[,-(2:6), drop = FALSE]

            hess[,1] <- rowSums(hess[,c(1,4,6)])
            hess <- hess[,-(2:6), drop = FALSE]
        }
    }

    colnames(grad) <- colnames(hess) <- param.names

    n.fixed <- length(fixed.param)
    if (n.fixed > 0){
        idx <- which(param.names %in% fixed.param)
        grad <- grad[, -idx, drop = FALSE]
        hess <- hess[, -idx, drop = FALSE]
    }

    if (any(is.na(grad)))
        return(list(var.score = NA, hessian = NA, gradient = NA))

    var.score <- var(grad) * n.obs
    hessian <- var(hess, na.rm = TRUE) * n.obs * n.pairs
    gradient <- as.double(colSums(grad))
    names(gradient) <- colnames(hessian)

    return(list(var.score = var.score, hessian = hessian, gradient = gradient))
}

.schlatherstderr <- function(par, data, dist, cov.mod, loc.dsgn.mat, scale.dsgn.mat,
                             shape.dsgn.mat, temp.dsgn.mat.loc, temp.dsgn.mat.scale,
                             temp.dsgn.mat.shape, use.temp.cov, fit.marge,
                             fixed.param, param.names, weights){

    ##data is a matrix with each column corresponds to one location
    ##distVec is the a matrix giving the "distance vector" for each pair
    ##(1 row = 1 station)
    n.site <- ncol(data)
    n.obs <- nrow(data)
    n.pairs <- n.site * (n.site - 1) / 2
    n.param <- length(param.names)

    if (is.null(weights))
        weights <- rep(1, n.pairs)

    nugget <- par["nugget"]
    range <- par["range"]
    smooth <- par["smooth"]

    if (cov.mod == 5)
        ##i.e. Generalized Cauchy
        smooth2 <- par["smooth2"]

    else
        ##it won't be used anyway...
        smooth2 <- 0

    if (fit.marge){

        n.loccoeff <- ncol(loc.dsgn.mat)
        n.scalecoeff <- ncol(scale.dsgn.mat)
        n.shapecoeff <- ncol(shape.dsgn.mat)

        loc.idx <- which(substr(names(par), 1, 3) == "loc")
        scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
        shape.idx <- which(substr(names(par), 1, 5) == "shape")

        loc.param <- par[loc.idx]
        scale.param <- par[scale.idx]
        shape.param <- par[shape.idx]

        n.temploccoeff <- ifelse(use.temp.cov[1], ncol(temp.dsgn.mat.loc), 0)
        n.tempscalecoeff <- ifelse(use.temp.cov[2], ncol(temp.dsgn.mat.scale), 0)
        n.tempshapecoeff <- ifelse(use.temp.cov[3], ncol(temp.dsgn.mat.scale), 0)

        temp.loc.idx <- which(substr(names(par), 1, 12) == "tempCoeffLoc")
        temp.scale.idx <- which(substr(names(par), 1, 14) == "tempCoeffScale")
        temp.shape.idx <- which(substr(names(par), 1, 14) == "tempCoeffShape")

        temp.loc.param <- ifelse(length(temp.loc.idx) > 0, par[temp.loc.idx], 0)
        temp.scale.param <- ifelse(length(temp.scale.idx) > 0, par[temp.scale.idx], 0)
        temp.shape.param <- ifelse(length(temp.shape.idx) > 0, par[temp.shape.idx], 0)
    }

    else {
        n.loccoeff <- n.scalecoeff <- n.shapecoeff <- loc.param <- scale.param <- shape.param <- 1
        temp.loc.param <- temp.scale.param <- temp.shape.param <- n.temploccoeff <-
            n.tempscalecoeff <- n.tempshapecoeff <- 0
    }

    std.err <- .C("schlatherstderr", as.integer(cov.mod), as.double(data), as.double(dist),
                  as.integer(n.site), as.integer(n.obs), as.double(loc.dsgn.mat),
                  as.integer(n.loccoeff), as.double(scale.dsgn.mat), as.integer(n.scalecoeff),
                  as.double(shape.dsgn.mat), as.integer(n.shapecoeff), as.double(temp.dsgn.mat.loc),
                  as.integer(n.temploccoeff), as.double(temp.dsgn.mat.scale),
                  as.integer(n.tempscalecoeff), as.double(temp.dsgn.mat.shape),
                  as.integer(n.tempshapecoeff), as.double(loc.param), as.double(scale.param),
                  as.double(shape.param), as.double(temp.loc.param), as.double(temp.scale.param),
                  as.double(temp.shape.param), as.double(nugget), as.double(range),
                  as.double(smooth), as.double(smooth2), as.integer(fit.marge),
                  as.integer(use.temp.cov), as.double(weights), hess = double(n.obs * n.param * n.pairs),
                  grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes", NAOK = TRUE)

    grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
    hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)
    colnames(grad) <- colnames(hess) <- param.names

    n.fixed <- length(fixed.param)
    if (n.fixed > 0){
        idx <- which(param.names %in% fixed.param)
        grad <- grad[, -idx, drop = FALSE]
        hess <- hess[, -idx, drop = FALSE]
    }

    if (any(is.na(grad)))
        return(list(var.score = NA, hessian = NA, gradient = NA))

    var.score <- var(grad) * n.obs
    hessian <- var(hess, na.rm = TRUE) * n.obs * n.pairs
    gradient <- as.double(colSums(grad))
    names(gradient) <- colnames(hessian)

    return(list(var.score = var.score, hessian = hessian, gradient = gradient))
}

.schlatherindstderr <- function(par, data, dist, cov.mod, loc.dsgn.mat, scale.dsgn.mat,
                                shape.dsgn.mat, temp.dsgn.mat.loc, temp.dsgn.mat.scale,
                                temp.dsgn.mat.shape, use.temp.cov, fit.marge,
                                fixed.param, param.names, weights){

    ##data is a matrix with each column corresponds to one location
    ##distVec is the a matrix giving the "distance vector" for each pair
    ##(1 row = 1 station)
    n.site <- ncol(data)
    n.obs <- nrow(data)
    n.pairs <- n.site * (n.site - 1) / 2
    n.param <- length(param.names)

    if (is.null(weights))
        weights <- rep(1, n.pairs)

    alpha <- par["alpha"]
    nugget <- par["nugget"]
    range <- par["range"]
    smooth <- par["smooth"]

    if (cov.mod == 5)
        ##i.e. Generalized cauchy
        smooth2 <- par["smooth2"]

    else
        ##It won't be used anyway
        smooth2 <- 0

    if (fit.marge){
        n.loccoeff <- ncol(loc.dsgn.mat)
        n.scalecoeff <- ncol(scale.dsgn.mat)
        n.shapecoeff <- ncol(shape.dsgn.mat)

        loc.idx <- which(substr(names(par), 1, 3) == "loc")
        scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
        shape.idx <- which(substr(names(par), 1, 5) == "shape")

        loc.param <- par[loc.idx]
        scale.param <- par[scale.idx]
        shape.param <- par[shape.idx]

        n.temploccoeff <- ifelse(use.temp.cov[1], ncol(temp.dsgn.mat.loc), 0)
        n.tempscalecoeff <- ifelse(use.temp.cov[2], ncol(temp.dsgn.mat.scale), 0)
        n.tempshapecoeff <- ifelse(use.temp.cov[3], ncol(temp.dsgn.mat.scale), 0)

        temp.loc.idx <- which(substr(names(par), 1, 12) == "tempCoeffLoc")
        temp.scale.idx <- which(substr(names(par), 1, 14) == "tempCoeffScale")
        temp.shape.idx <- which(substr(names(par), 1, 14) == "tempCoeffShape")

        temp.loc.param <- ifelse(length(temp.loc.idx) > 0, par[temp.loc.idx], 0)
        temp.scale.param <- ifelse(length(temp.scale.idx) > 0, par[temp.scale.idx], 0)
        temp.shape.param <- ifelse(length(temp.shape.idx) > 0, par[temp.shape.idx], 0)
    }

    else {
        n.loccoeff <- n.scalecoeff <- n.shapecoeff <- loc.param <- scale.param <- shape.param <- 1
        temp.loc.param <- temp.scale.param <- temp.shape.param <- n.temploccoeff <-
            n.tempscalecoeff <- n.tempshapecoeff <- 0
    }

    std.err <- .C("schlatherindstderr", as.integer(cov.mod), as.double(data), as.double(dist),
                  as.integer(n.site), as.integer(n.obs), as.double(loc.dsgn.mat),
                  as.integer(n.loccoeff), as.double(scale.dsgn.mat), as.integer(n.scalecoeff),
                  as.double(shape.dsgn.mat), as.integer(n.shapecoeff),
                  as.double(temp.dsgn.mat.loc), as.integer(n.temploccoeff),
                  as.double(temp.dsgn.mat.scale), as.integer(n.tempscalecoeff),
                  as.double(temp.dsgn.mat.shape), as.integer(n.tempshapecoeff),
                  as.double(loc.param), as.double(scale.param), as.double(shape.param),
                  as.double(temp.loc.param), as.double(temp.scale.param),
                  as.double(temp.shape.param), as.double(alpha), as.double(nugget),
                  as.double(range), as.double(smooth), as.double(smooth2), as.integer(fit.marge),
                  as.integer(use.temp.cov), as.double(weights), hess = double(n.obs * n.param * n.pairs),
                  grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes", NAOK = TRUE)

    grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
    hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)
    colnames(grad) <- colnames(hess) <- param.names

    n.fixed <- length(fixed.param)
    if (n.fixed > 0){
        idx <- which(param.names %in% fixed.param)
        grad <- grad[, -idx, drop = FALSE]
        hess <- hess[, -idx, drop = FALSE]
    }

    if (any(is.na(grad)))
        return(list(var.score = NA, hessian = NA, gradient = NA))

    var.score <- var(grad) * n.obs
    hessian <- var(hess, na.rm = TRUE) * n.obs * n.pairs
    gradient <- as.double(colSums(grad))
    names(gradient) <- colnames(hessian)

    return(list(var.score = var.score, hessian = hessian, gradient = gradient))
}

.geomgaussstderr <- function(par, data, dist, cov.mod, loc.dsgn.mat, scale.dsgn.mat,
                             shape.dsgn.mat, temp.dsgn.mat.loc, temp.dsgn.mat.scale,
                             temp.dsgn.mat.shape, use.temp.cov, fit.marge,
                             fixed.param, param.names, weights){

    ##data is a matrix with each column corresponds to one location
    ##distVec is the a matrix giving the "distance vector" for each pair
    ##(1 row = 1 station)
    n.site <- ncol(data)
    n.obs <- nrow(data)
    n.pairs <- n.site * (n.site - 1) / 2
    n.param <- length(param.names)

    if (is.null(weights))
        weights <- rep(1, n.pairs)

    sigma2 <- par["sigma2"]
    nugget <- par["nugget"]
    range <- par["range"]
    smooth <- par["smooth"]

    if (cov.mod == 5)
        ##i.e. Generalized cauchy
        smooth2 <- par["smooth2"]

    else
        ##it won't be used anyway
        smooth2 <- 0

    if (fit.marge){

        n.loccoeff <- ncol(loc.dsgn.mat)
        n.scalecoeff <- ncol(scale.dsgn.mat)
        n.shapecoeff <- ncol(shape.dsgn.mat)

        loc.idx <- which(substr(names(par), 1, 3) == "loc")
        scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
        shape.idx <- which(substr(names(par), 1, 5) == "shape")

        loc.param <- par[loc.idx]
        scale.param <- par[scale.idx]
        shape.param <- par[shape.idx]

        n.temploccoeff <- ifelse(use.temp.cov[1], ncol(temp.dsgn.mat.loc), 0)
        n.tempscalecoeff <- ifelse(use.temp.cov[2], ncol(temp.dsgn.mat.scale), 0)
        n.tempshapecoeff <- ifelse(use.temp.cov[3], ncol(temp.dsgn.mat.scale), 0)

        temp.loc.idx <- which(substr(names(par), 1, 12) == "tempCoeffLoc")
        temp.scale.idx <- which(substr(names(par), 1, 14) == "tempCoeffScale")
        temp.shape.idx <- which(substr(names(par), 1, 14) == "tempCoeffShape")

        temp.loc.param <- ifelse(length(temp.loc.idx) > 0, par[temp.loc.idx], 0)
        temp.scale.param <- ifelse(length(temp.scale.idx) > 0, par[temp.scale.idx], 0)
        temp.shape.param <- ifelse(length(temp.shape.idx) > 0, par[temp.shape.idx], 0)
    }

    else {
        n.loccoeff <- n.scalecoeff <- n.shapecoeff <- loc.param <- scale.param <- shape.param <- 1
        temp.loc.param <- temp.scale.param <- temp.shape.param <- n.temploccoeff <-
            n.tempscalecoeff <- n.tempshapecoeff <- 0
    }

    std.err <- .C("geomgaussstderr", as.integer(cov.mod), as.double(data), as.double(dist),
                  as.integer(n.site), as.integer(n.obs), as.double(loc.dsgn.mat),
                  as.integer(n.loccoeff), as.double(scale.dsgn.mat), as.integer(n.scalecoeff),
                  as.double(shape.dsgn.mat), as.integer(n.shapecoeff),
                  as.double(temp.dsgn.mat.loc), as.integer(n.temploccoeff),
                  as.double(temp.dsgn.mat.scale), as.integer(n.tempscalecoeff),
                  as.double(temp.dsgn.mat.shape), as.integer(n.tempshapecoeff),
                  as.double(loc.param), as.double(scale.param), as.double(shape.param),
                  as.double(temp.loc.param), as.double(temp.scale.param),
                  as.double(temp.shape.param), as.double(sigma2), as.double(nugget),
                  as.double(range), as.double(smooth), as.double(smooth2), as.integer(fit.marge),
                  as.integer(use.temp.cov), as.double(weights), hess = double(n.obs * n.param * n.pairs),
                  grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes", NAOK = TRUE)

    grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
    hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)
    colnames(grad) <- colnames(hess) <- param.names

    n.fixed <- length(fixed.param)
    if (n.fixed > 0){
        idx <- which(param.names %in% fixed.param)
        grad <- grad[, -idx, drop = FALSE]
        hess <- hess[, -idx, drop = FALSE]
    }

    if (any(is.na(grad)))
        return(list(var.score = NA, hessian = NA, gradient = NA))

    var.score <- var(grad) * n.obs
    hessian <- var(hess, na.rm = TRUE) * n.obs * n.pairs
    gradient <- as.double(colSums(grad))
    names(gradient) <- colnames(hessian)

    return(list(var.score = var.score, hessian = hessian, gradient = gradient))
}


.brownresnickstderr <- function(par, data, dist, loc.dsgn.mat, scale.dsgn.mat,
                                shape.dsgn.mat, temp.dsgn.mat.loc, temp.dsgn.mat.scale,
                                temp.dsgn.mat.shape, use.temp.cov, fit.marge,
                                fixed.param, param.names, weights){

    ##data is a matrix with each column corresponds to one location
    ##distVec is the a matrix giving the "distance vector" for each pair
    ##(1 row = 1 station)
    n.site <- ncol(data)
    n.obs <- nrow(data)
    n.pairs <- n.site * (n.site - 1) / 2
    n.param <- length(param.names)

    if (is.null(weights))
        weights <- rep(1, n.pairs)

    range <- par["range"]
    smooth <- par["smooth"]

    if (fit.marge){

        n.loccoeff <- ncol(loc.dsgn.mat)
        n.scalecoeff <- ncol(scale.dsgn.mat)
        n.shapecoeff <- ncol(shape.dsgn.mat)

        loc.idx <- which(substr(names(par), 1, 3) == "loc")
        scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
        shape.idx <- which(substr(names(par), 1, 5) == "shape")

        loc.param <- par[loc.idx]
        scale.param <- par[scale.idx]
        shape.param <- par[shape.idx]

        n.temploccoeff <- ifelse(use.temp.cov[1], ncol(temp.dsgn.mat.loc), 0)
        n.tempscalecoeff <- ifelse(use.temp.cov[2], ncol(temp.dsgn.mat.scale), 0)
        n.tempshapecoeff <- ifelse(use.temp.cov[3], ncol(temp.dsgn.mat.scale), 0)

        temp.loc.idx <- which(substr(names(par), 1, 12) == "tempCoeffLoc")
        temp.scale.idx <- which(substr(names(par), 1, 14) == "tempCoeffScale")
        temp.shape.idx <- which(substr(names(par), 1, 14) == "tempCoeffShape")

        temp.loc.param <- ifelse(length(temp.loc.idx) > 0, par[temp.loc.idx], 0)
        temp.scale.param <- ifelse(length(temp.scale.idx) > 0, par[temp.scale.idx], 0)
        temp.shape.param <- ifelse(length(temp.shape.idx) > 0, par[temp.shape.idx], 0)
    }

    else {
        n.loccoeff <- n.scalecoeff <- n.shapecoeff <- loc.param <- scale.param <- shape.param <- 1
        temp.loc.param <- temp.scale.param <- temp.shape.param <- n.temploccoeff <-
            n.tempscalecoeff <- n.tempshapecoeff <- 0
    }

    std.err <- .C("brownresnickstderr", as.double(data), as.double(dist), as.integer(n.site),
                  as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
                  as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                  as.integer(n.shapecoeff), as.double(temp.dsgn.mat.loc),
                  as.integer(n.temploccoeff), as.double(temp.dsgn.mat.scale),
                  as.integer(n.tempscalecoeff), as.double(temp.dsgn.mat.shape),
                  as.integer(n.tempshapecoeff), as.double(loc.param), as.double(scale.param),
                  as.double(shape.param), as.double(temp.loc.param), as.double(temp.scale.param),
                  as.double(temp.shape.param), as.double(range), as.double(smooth),
                  as.integer(fit.marge), as.integer(use.temp.cov), as.double(weights),
                  hess = double(n.obs * n.param * n.pairs), grad = double(n.obs * n.param),
                  PACKAGE = "SpatialExtremes", NAOK = TRUE)

    grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
    hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)
    colnames(grad) <- colnames(hess) <- param.names

    n.fixed <- length(fixed.param)
    if (n.fixed > 0){
        idx <- which(param.names %in% fixed.param)
        grad <- grad[, -idx, drop = FALSE]
        hess <- hess[, -idx, drop = FALSE]
    }

    if (any(is.na(grad)))
        return(list(var.score = NA, hessian = NA, gradient = NA))

    var.score <- var(grad) * n.obs
    hessian <- var(hess, na.rm = TRUE) * n.obs * n.pairs
    gradient <- as.double(colSums(grad))
    names(gradient) <- colnames(hessian)

    return(list(var.score = var.score, hessian = hessian, gradient = gradient))
}

.spatgevstderr <- function(par, data, loc.dsgn.mat, scale.dsgn.mat, shape.dsgn.mat,
                           temp.dsgn.mat.loc, temp.dsgn.mat.scale, temp.dsgn.mat.shape,
                           use.temp.cov, fixed.param, param.names){

    ##data is a matrix with each column corresponds to one location
    n.site <- ncol(data)
    n.obs <- nrow(data)
    n.param <- length(param.names)

    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)

    loc.idx <- which(substr(names(par), 1, 3) == "loc")
    scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
    shape.idx <- which(substr(names(par), 1, 5) == "shape")

    loc.param <- par[loc.idx]
    scale.param <- par[scale.idx]
    shape.param <- par[shape.idx]

    n.temploccoeff <- ifelse(use.temp.cov[1], ncol(temp.dsgn.mat.loc), 0)
    n.tempscalecoeff <- ifelse(use.temp.cov[2], ncol(temp.dsgn.mat.scale), 0)
    n.tempshapecoeff <- ifelse(use.temp.cov[3], ncol(temp.dsgn.mat.scale), 0)

    temp.loc.idx <- which(substr(names(par), 1, 12) == "tempCoeffLoc")
    temp.scale.idx <- which(substr(names(par), 1, 14) == "tempCoeffScale")
    temp.shape.idx <- which(substr(names(par), 1, 14) == "tempCoeffShape")

    temp.loc.param <- ifelse(length(temp.loc.idx) > 0, par[temp.loc.idx], 0)
    temp.scale.param <- ifelse(length(temp.scale.idx) > 0, par[temp.scale.idx], 0)
    temp.shape.param <- ifelse(length(temp.shape.idx) > 0, par[temp.shape.idx], 0)

    std.err <- .C("spatgevstderr", as.double(data), as.integer(n.site), as.integer(n.obs),
                  as.double(loc.dsgn.mat), as.integer(n.loccoeff), as.double(scale.dsgn.mat),
                  as.integer(n.scalecoeff), as.double(shape.dsgn.mat), as.integer(n.shapecoeff),
                  as.double(temp.dsgn.mat.loc), as.integer(n.temploccoeff),
                  as.double(temp.dsgn.mat.scale), as.integer(n.tempscalecoeff),
                  as.double(temp.dsgn.mat.shape), as.integer(n.tempshapecoeff),
                  as.double(loc.param), as.double(scale.param), as.double(shape.param),
                  as.double(temp.loc.param), as.double(temp.scale.param),
                  as.double(temp.shape.param), as.integer(use.temp.cov),
                  hess = double(n.obs * n.param * n.site), grad = double(n.obs * n.param),
                  PACKAGE = "SpatialExtremes", NAOK = TRUE)

    grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
    hess <- matrix(std.err$hess, nrow = n.obs * n.site, ncol = n.param)
    colnames(grad) <- colnames(hess) <- param.names

    n.fixed <- length(fixed.param)
    if (n.fixed > 0){
        idx <- which(param.names %in% fixed.param)
        hess <- hess[, -idx, drop = FALSE]
        grad <- grad[, -idx, drop = FALSE]
    }

    if (any(is.na(grad)))
        return(list(gradient = NA, var.score = NA, hessian = NA))

    var.score <- var(grad) * n.obs
    hessian <- var(hess, na.rm = TRUE) * n.obs * n.site

    gradient <- as.double(colSums(grad))
    names(gradient) <- colnames(hessian)

    return(list(var.score = var.score, hessian = hessian, gradient = gradient))
}

.extremaltstderr <- function(par, data, dist, cov.mod, loc.dsgn.mat, scale.dsgn.mat,
                             shape.dsgn.mat, temp.dsgn.mat.loc, temp.dsgn.mat.scale,
                             temp.dsgn.mat.shape, use.temp.cov, fit.marge,
                             fixed.param, param.names, weights){

    ##data is a matrix with each column corresponds to one location
    ##distVec is the a matrix giving the "distance vector" for each pair
    ##(1 row = 1 station)
    n.site <- ncol(data)
    n.obs <- nrow(data)
    n.pairs <- n.site * (n.site - 1) / 2
    n.param <- length(param.names)

    if (is.null(weights))
        weights <- rep(1, n.pairs)

    nugget <- par["nugget"]
    range <- par["range"]
    smooth <- par["smooth"]
    DoF <- par["DoF"]

    if (cov.mod == 5)
        ##i.e. Generalized Cauchy
        smooth2 <- par["smooth2"]

    else
        ##it won't be used anyway...
        smooth2 <- 0

    if (fit.marge){

        n.loccoeff <- ncol(loc.dsgn.mat)
        n.scalecoeff <- ncol(scale.dsgn.mat)
        n.shapecoeff <- ncol(shape.dsgn.mat)

        loc.idx <- which(substr(names(par), 1, 3) == "loc")
        scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
        shape.idx <- which(substr(names(par), 1, 5) == "shape")

        loc.param <- par[loc.idx]
        scale.param <- par[scale.idx]
        shape.param <- par[shape.idx]

        n.temploccoeff <- ifelse(use.temp.cov[1], ncol(temp.dsgn.mat.loc), 0)
        n.tempscalecoeff <- ifelse(use.temp.cov[2], ncol(temp.dsgn.mat.scale), 0)
        n.tempshapecoeff <- ifelse(use.temp.cov[3], ncol(temp.dsgn.mat.scale), 0)

        temp.loc.idx <- which(substr(names(par), 1, 12) == "tempCoeffLoc")
        temp.scale.idx <- which(substr(names(par), 1, 14) == "tempCoeffScale")
        temp.shape.idx <- which(substr(names(par), 1, 14) == "tempCoeffShape")

        temp.loc.param <- ifelse(length(temp.loc.idx) > 0, par[temp.loc.idx], 0)
        temp.scale.param <- ifelse(length(temp.scale.idx) > 0, par[temp.scale.idx], 0)
        temp.shape.param <- ifelse(length(temp.shape.idx) > 0, par[temp.shape.idx], 0)
    }

    else {
        n.loccoeff <- n.scalecoeff <- n.shapecoeff <- loc.param <- scale.param <- shape.param <- 1
        temp.loc.param <- temp.scale.param <- temp.shape.param <- n.temploccoeff <-
            n.tempscalecoeff <- n.tempshapecoeff <- 0
    }

    std.err <- .C("extremaltstderr", as.integer(cov.mod), as.double(data), as.double(dist),
                  as.integer(n.site), as.integer(n.obs), as.double(loc.dsgn.mat),
                  as.integer(n.loccoeff), as.double(scale.dsgn.mat), as.integer(n.scalecoeff),
                  as.double(shape.dsgn.mat), as.integer(n.shapecoeff), as.double(temp.dsgn.mat.loc),
                  as.integer(n.temploccoeff), as.double(temp.dsgn.mat.scale),
                  as.integer(n.tempscalecoeff), as.double(temp.dsgn.mat.shape),
                  as.integer(n.tempshapecoeff), as.double(loc.param), as.double(scale.param),
                  as.double(shape.param), as.double(temp.loc.param), as.double(temp.scale.param),
                  as.double(temp.shape.param), as.double(nugget), as.double(range),
                  as.double(smooth), as.double(smooth2), as.double(DoF), as.integer(fit.marge),
                  as.integer(use.temp.cov), as.double(weights), hess = double(n.obs * n.param * n.pairs),
                  grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes", NAOK = TRUE)

    grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
    hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)
    colnames(grad) <- colnames(hess) <- param.names

    n.fixed <- length(fixed.param)
    if (n.fixed > 0){
        idx <- which(param.names %in% fixed.param)
        grad <- grad[, -idx, drop = FALSE]
        hess <- hess[, -idx, drop = FALSE]
    }

    if (any(is.na(grad)))
        return(list(var.score = NA, hessian = NA, gradient = NA))

    var.score <- var(grad) * n.obs
    hessian <- var(hess, na.rm = TRUE) * n.obs * n.pairs
    gradient <- as.double(colSums(grad))
    names(gradient) <- colnames(hessian)

    return(list(var.score = var.score, hessian = hessian, gradient = gradient))
}
