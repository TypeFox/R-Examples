latent <- function(data, coord, cov.mod = "powexp", loc.form, scale.form,
                   shape.form, marg.cov = NULL, hyper, prop, start, n = 5000,
                   thin = 1, burn.in = 0){

    if (!(all(cov.mod %in% c("whitmat","cauchy","powexp","bessel"))))
        stop("''cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp', 'bessel'")

    if (!(length(cov.mod) %in% c(1,3)))
        stop("''cov.mod'' must be of length 1 or 3")

    if (!is.list(hyper))
        stop("''hyper'' must be a list")

    if (is.null(hyper$sills) || is.null(hyper$ranges) || is.null(hyper$smooths) ||
        is.null(hyper$betaMeans) || is.null(hyper$betaIcov))
        stop("''hyper'' must be a named list with named components 'sills', 'ranges',
 'smooths', 'betaMeans' and 'betaIcov'")

    for (i in 1:length(hyper))
        if (!all(names(hyper[[i]]) %in% c("loc", "scale", "shape")) ||
            length(hyper[[i]]) != 3)
            stop("Each component of the list ''hyper'' must a be a named list with names
 'loc', 'scale' and 'shape'")

    if (is.null(prop$gev) || is.null(prop$ranges) || is.null(prop$smooths))
        stop("''prop'' must be a named list with named components 'gev', 'ranges', 'smooths'")

    for (i in 1:3)
        if (length(prop[[i]]) != 3)
            stop("Each component of the named list 'prop' must be a vector of length 3 ---
 one value for each margin")

    if (!is.list(start))
        stop("''start'' must be a list")

    if (is.null(start$sills) || is.null(start$ranges) || is.null(start$smooths) ||
        is.null(start$beta))
        stop("''start'' must be a named list with named components 'sills', 'ranges', 'smooths'
 and 'beta'")

    if (is.null(start$beta$loc) || is.null(start$beta$scale) || is.null(start$beta$shape))
        stop("''start$beta'' must be a named list with named components 'loc', 'scale' and 'shape'")

    for (i in 1:3)
        if ((length(start$sills) != 3) || (length(start$ranges) != 3) ||
            (length(start$smooths) != 3))
            stop("''start$sills'', ''start$ranges'' and ''start$smooths'' must be numeric vectors of length 3")

    if (length(cov.mod) == 1)
        cov.mod <- rep(cov.mod, 3)

    cov.mod.num <- rep(NA, 3)

    for (i in 1:3){
        if (cov.mod[i] == "whitmat")
            cov.mod.num[i] <- 1
        if (cov.mod[i] == "cauchy")
            cov.mod.num[i] <- 2
        if (cov.mod[i] == "powexp")
            cov.mod.num[i] <- 3
        if (cov.mod[i] == "bessel")
            cov.mod.num[i] <- 4
    }

    n.site <- ncol(data)
    n.obs <- nrow(data)
    dist.dim <- ncol(coord)
    distMat <- t(as.matrix(dist(coord, diag = TRUE)))
    distMat <- distMat[lower.tri(distMat, diag = TRUE)]

    ##Only the upper diagonal elements have to be stored for the
    ##hyperparameters --- setting the lower diagonal elements to 0
    hyper$betaIcov$loc[lower.tri(hyper$betaIcov$loc)] <- 0
    hyper$betaIcov$scale[lower.tri(hyper$betaIcov$scale)] <- 0
    hyper$betaIcov$shape[lower.tri(hyper$betaIcov$shape)] <- 0

    ##With our notation, formula must be of the form y ~ xxxx
    loc.form <- update(loc.form, y ~ .)
    scale.form <- update(scale.form, y ~ .)
    shape.form <- update(shape.form, y ~ .)

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

    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)
    n.beta <- c(n.loccoeff, n.scalecoeff, n.shapecoeff)

    for (i in 1:3)
        if ((length(hyper$betaMeans[[i]]) != n.beta[i]) ||
            (any(dim(hyper$betaIcov[[i]]) != c(n.beta[i], n.beta[i]))))
            stop("The hyper parameters for the regression parameters doesn't match")

    for (i in 1:3)
        if (length(start$beta[[i]]) != n.beta[i])
            stop("The starting values for the regression parameters doesn't match")

    GEVmles <- t(apply(data, 2, gevmle))

    chain.loc <- double(n * (n.loccoeff + 3 + n.site))
    chain.scale <- double(n * (n.scalecoeff + 3 + n.site))
    chain.shape <- double(n * (n.shapecoeff + 3 + n.site))

    temp <- .C("latentgev", as.integer(n), as.double(data), as.integer(n.site),
               as.integer(n.obs), as.integer(cov.mod.num), as.integer(dist.dim),
               as.double(distMat), as.double(c(loc.dsgn.mat, scale.dsgn.mat, shape.dsgn.mat)),
               as.integer(n.beta), as.double(unlist(start$beta)), as.double(start$sills),
               as.double(start$ranges), as.double(start$smooths),
               as.double(GEVmles), as.double(unlist(hyper$sills)),
               as.double(unlist(hyper$ranges)), as.double(unlist(hyper$smooths)),
               as.double(unlist(hyper$betaMeans)), as.double(unlist(hyper$betaIcov)),
               as.double(prop$gev), as.double(prop$ranges), as.double(prop$smooths),
               chain.loc = chain.loc, chain.scale = chain.scale,  chain.shape = chain.shape,
               acc.rates = double(9), ext.rates = double(9), as.integer(thin),
               burn.in = as.integer(burn.in), PACKAGE = "SpatialExtremes", NAOK = TRUE)

    chain.loc <- matrix(temp$chain.loc, n, byrow = TRUE)
    chain.scale <- matrix(temp$chain.scale, n, byrow = TRUE)
    chain.shape <- matrix(temp$chain.shape, n, byrow = TRUE)
    acc.rates <- temp$acc.rates
    ext.rates <- temp$ext.rates
    names(acc.rates) <- names(ext.rates) <-
        c("gev:loc", "gev:scale", "gev:shape", "range:loc", "range:scale", "range:shape",
          "smooth:loc", "smooth:scale", "smooth:shape")

    acc.rates[1:3] <- acc.rates[1:3] / n.site
    ext.rates[1:3] <- ext.rates[1:3] / n.site

    colnames(chain.loc) <- c(paste("lm", 1:n.loccoeff,sep=""),
                             "sill", "range", "smooth",
                             paste("loc", 1:n.site,sep=""))
    colnames(chain.scale) <- c(paste("lm", 1:n.scalecoeff,sep=""),
                               "sill", "range", "smooth",
                               paste("scale", 1:n.site,sep=""))
    colnames(chain.shape) <- c(paste("lm", 1:n.shapecoeff,sep=""),
                               "sill", "range", "smooth",
                               paste("shape", 1:n.site,sep=""))

    mcmc <- list(chain.loc = chain.loc, chain.scale = chain.scale,
                 chain.shape = chain.shape, loc.dsgn.mat = loc.dsgn.mat,
                 scale.dsgn.mat = scale.dsgn.mat, shape.dsgn.mat = shape.dsgn.mat,
                 acc.rates = rbind(acc.rates = acc.rates, ext.rates = ext.rates),
                 hyper = hyper, cov.mod = cov.mod, burn.in = burn.in, thin = thin,
                 data = data, coord = coord, marg.cov = marg.cov, loc.form = loc.form,
                 scale.form = scale.form, shape.form = shape.form)
    class(mcmc) <- "latent"
    dummy <- DIC(mcmc)
    mcmc <- c(mcmc, list(eNoP = dummy["eNoP"], DIC = dummy["DIC"],
                         Dbar = dummy["Dbar"]))
    class(mcmc) <- "latent"
    return(mcmc)
}

DIC <- function(object, ..., fun = "mean"){
    if (class(object) != "latent")
        stop("'DIC' can only be used with objects of class 'latent'")

    chain.loc <- object$chain.loc
    chain.scale <- object$chain.scale
    chain.shape <- object$chain.shape

    loc.idx <- which(substr(colnames(chain.loc), 1, 3) == "loc")
    scale.idx <- which(substr(colnames(chain.scale), 1, 5) == "scale")
    shape.idx <- which(substr(colnames(chain.shape), 1, 5) == "shape")

    chain.loc <- chain.loc[,loc.idx]
    chain.scale <- chain.scale[,scale.idx]
    chain.shape <- chain.shape[,shape.idx]

    n.chain <- nrow(chain.loc)
    n.site <- ncol(object$data)
    n.obs <- nrow(object$data)
    post.loc <- apply(chain.loc, 2, fun)
    post.scale <- apply(chain.scale, 2, fun)
    post.shape <- apply(chain.shape, 2, fun)

    temp <- .C("DIC", as.integer(n.chain), as.integer(n.site), as.integer(n.obs),
               as.double(object$data), as.double(chain.loc),
               as.double(chain.scale), as.double(chain.shape),
               as.double(post.loc), as.double(post.scale), as.double(post.shape),
               dic = double(1), effNpar = double(1), dbar = double(1),
               PACKAGE = "SpatialExtremes", NAOK = TRUE)

    return(c(DIC = temp$dic, eNoP = temp$effNpar, Dbar = temp$dbar))
}


