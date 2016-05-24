#require(MASS)
#require(lavaan)
#require(numDeriv)

## Mirror the lower or upper diagonal of a square matrix
mirror.tri <- function(mat, upper.on.lower = FALSE) {
    if (upper.on.lower) {
        mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    } else {
        mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
    }
    mat
}

## Factor analysis using principal components as the extraction method
factanal.prcomp <- function (covmat, nfactors = NULL) {
    prc <- princomp(covmat = covmat, cor = TRUE)
    eig <- prc$sdev^2
    if (is.null(nfactors)) nfactors <- sum(eig >= 1) 
    loads <- prc$loadings[, 1:nfactors]
    loads[, 1:nfactors] %*% diag(prc$sdev[1:nfactors])
}

## Get a covariance matrix from a correlation matrix and std deviations
cor2cov <- function(cormat, sds) {
    outer(sds, sds, "*") * cormat
}

## Check if a square matrix is Hermitian
is.hermitian <- function(covmat) {
    t.covmat <- t(covmat)
    if (all(dim(covmat) == dim(t.covmat)) && all(covmat == t(covmat))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

## Check if all covariances are within bounds
is.within.bounds <- function(covmat) {
    if (!is.hermitian(covmat)) stop("Matrix is not Hermitian")
    sds <- sqrt(diag(covmat))
    max.cov <- outer(sds, sds, "*")
    differences <- round(max.cov - covmat, 3)
    res <- all(differences >= 0)
    if (!res) {
      outside.bounds <- which(differences < 0, arr.ind = TRUE)
      varnames <- rownames(differences)
      cat("Out of bounds covariances:\n\n")
      if (!is.null(varnames)) {
          for (i in 1:(NROW(outside.bounds)/2)) {
              rn <- varnames[outside.bounds[i, 1]]
              cn <- varnames[outside.bounds[i, 2]]
              cat(sprintf("%s and %s: max - %f, obs - %f", rn,
                          cn, max.cov[rn, cn], covmat[rn, cn]), "\n") 
          }
      } else {
          for (i in 1:(NROW(outside.bounds)/2)) {
              rn <- outside.bounds[i, 1]
              cn <- outside.bounds[i, 2]
              cat(sprintf("row %s and column %s: max - %f, obs - %f", rn,
                          cn, max.cov[rn, cn], covmat[rn, cn]), "\n") 
          }
      }
      cat("\n")
    }
    res
}

## Check if a covariance matrix is positive-definite
## We could try to invert it, but instead we opted to do a bunch of
## testing for learning purposes.
is.pd <- function(covmat) {
    ## Check if matrix is Hermitian (lower and upper diag are equal)
    is.ht <- is.hermitian(covmat)
    if (!is.ht) {
        warning("Matrix isn't Hermitian")
        return(FALSE)
    }
    ## Check if all covariances are within bounds
    is.wd <- is.within.bounds(covmat)
    if (!is.wd) {
        warning("Covariances aren't within bounds")
        return(FALSE)
    }
    ## Check if the matrix is invertible (this whole funcion could be
    ## replaced by a call to chol
    is.invertible <- try(solve(covmat), silent = TRUE)
    if (class(is.invertible) == "try-error") {
        warning("Matrix isn't invertible")
        return(FALSE)
    }
    ## Check if all eigenvalues are positive
    eigenvalues <- eigen(covmat)$values
    if (!all(eigenvalues > 0)) {
        warning("Matrix has negative eigenvalues")
        return(FALSE)
    }
    ## Check if the determinant is higher than 0
    determ <- prod(eigenvalues)
    if (determ <= 0) {
        warning("The determinant is equal or lower than 0")
        return(FALSE)
    }
    TRUE
}


## Creates a multivariate normal distribution from a covariance matrix.
rmnorm <- function(covmat, means = 0, sds = NULL, n) {
    if (is.null(sds)) 
        sds <- sqrt(diag(covmat))
    pattern <- factanal.prcomp(covmat = covmat, nfactors = NROW(covmat))
    nvar <- NROW(covmat)
    X <- sapply(seq_len(nvar), function(x) {rnorm(n)})
    X <- t(X)
    Z <- pattern %*% X
    Z <- t(Z)
    Z <- sweep(Z, 2, sds, "*")
    Z <- sweep(Z, 2, means, "+")
    Z
}

## Identify the type of a lavaanified model
model.type <- function(model) {
    if(all(model$op != "=~")) 
        type <- "pa"
    if(any(model$op == "=~") & any(model$op == "~")) 
        type <- "sem"
    if(any(model$op == "=~") & all(model$op != "~")) 
        type <- "cfa"
    type
}

## Build but dont fit a model using lavaan to use as a starting point.
lvnfy <- function(model, data, n = NULL, ...) {
    if (is.hermitian(data)) {
        out <- sem(model, sample.cov = data, sample.nobs = n, do.fit
                      = FALSE, ...)
    } else {
        out <- sem(model, data = data, do.fit = FALSE, ...)
    }
    lv.list <- inspect(out, "list")
    lv.list$ustart <-out@Fit@start
    lv.list <- as.list(lv.list[c(1:4,7:8)])
    vars <- variables(lv.list)
    var.names <- rownames(vars)
    latents <- var.names[!vars[,"lat"]]
    regressed <- !vars[,"lat"] & (vars[,"endo"] | vars[,"exo"])
    vars <- cbind(vars, "latenterise" = regressed)
    attr(lv.list, "vars") <- vars
    lv.list
}

## Get latent variable names
latents <- function(model) {
    vars <- attr(model, "vars")
    rownames(vars)[vars[,"lat"]]
}

## Get observed variables (even if they are not indicators)
observed <- function(model) {
    vars <- attr(model, "vars")
    rownames(vars)[!vars[,"lat"]]
}

## Get indicators and separate them by type
manifests <- function(model) {
    vars <- attr(model, "vars")
    rownames(vars)[vars[,"ind"]]
}

## Get endogenous variables
endogenous <- function(model) {
    vars <- attr(model, "vars")
    rownames(vars)[vars[,"endo"]]
}

## Get endogenous variables
exogenous <- function(model) {
    vars <- attr(model, "vars")
    rownames(vars)[vars[,"exo"]]
}

## Identity (help done)
II <- function(model) {
    II <- BETA(model)
    II[] <- 0
    diag(II) <- 1
    II
}

## Beta matrix (regression coefficients between endogenous variables)
BETA <- function(model) {
    vars <- attr(model, "vars")
    var.names <- rownames(vars)
    ## latents
    latents <- vars[,"lat"]
    obs.not.indicators <- !vars[,"lat"] & (vars[,"exo"] | vars[,"endo"])
    sides <- var.names[latents | obs.not.indicators]
    regressions <- lapply(model, `[`, model$op == "~" | (model$op ==
                                                         "=~" &
                                                         model$rhs %in%
                                                         var.names[latents]))
    sides <- unique(var.names[latents | obs.not.indicators])
    beta.mat <- matrix(0, length(sides), length(sides))
    dimnames(beta.mat) <- list(sides, sides)
    free.pars.pos <- NULL
    index.matrix <- beta.mat
    index.matrix[] <- seq_along(index.matrix)
    free.pars <- regressions$free[regressions$free > 0L]
    if (length(regressions) > 0) {
        for (i in seq_len(length(regressions$lhs))) {
            lhs <- regressions$lhs[i]
            rhs <- regressions$rhs[i]
            op <- regressions$op[i]
            start.value <- regressions$ustart[i] 
            is.free <- regressions$free[i] > 0L
            if (op == "~") {
                beta.mat[lhs, rhs] <- start.value
                if (is.free) {
                    free.pars.pos <- c(free.pars.pos, index.matrix[lhs,
                                       rhs])
                }
            } else{
                beta.mat[rhs, lhs] <- start.value
                if (is.free) {
                    free.pars.pos <- c(free.pars.pos, index.matrix[rhs,
                                       lhs])
                }
            }
        }
    }
    attr(beta.mat, "free.pars") <- free.pars 
    attr(beta.mat, "free.pars.pos") <- free.pars.pos
    beta.mat
}

## Psi matrix (variances and covariances of observed variables).
PSI <- function(model) {
    vars <- attr(model, "vars")
    var.names <- rownames(vars)
    latents <- vars[,"lat"]
    obs.not.indicators <- !vars[,"lat"] & (vars[,"exo"] | vars[,"endo"])
    sides <- var.names[latents | obs.not.indicators]
    psi.mat <- matrix(0, length(sides), length(sides))
    dimnames(psi.mat) <- list(sides, sides)
    var.covar <- lapply(model, `[`, model$op == "~~" & model$lhs %in%
                        sides & model$rhs %in% sides)
    #free.pars.pos <- numeric(length(var.covar$lhs) * 2 - NCOL(psi.mat))
    #free.pars <- numeric(length(free.pars.pos))
    free.pars.pos <- free.pars <- NULL
    index.matrix <- psi.mat
    index.matrix[] <- seq_along(index.matrix)
    if (length(var.covar) == 0) stop("Looks like you forgot to explicitly
                                     assign at least the disturbances in the model.")
    for (i in seq_len(length(var.covar$lhs))) {
        lhs <- var.covar$lhs[i]
        rhs <- var.covar$rhs[i]
        start.value <- var.covar$ustart[i]
        pos <- var.covar$free[i]
        is.free <- pos > 0L
        if (lhs == rhs) {
            psi.mat[lhs, rhs] <- start.value
            if (is.free) {
                free.pars.pos <- c(free.pars.pos, index.matrix[lhs, rhs])
                free.pars <- c(free.pars, pos)
            }
        } else {
            psi.mat[lhs, rhs] <- psi.mat[rhs, lhs] <- start.value
            if (is.free) {
                free.pars.pos <- c(free.pars.pos, index.matrix[lhs, rhs],
                                   index.matrix[rhs, lhs])
                free.pars <- c(free.pars, rep(pos, 2))
            }
        }
    }
    if (length(psi.mat) == 0L) {
        return(0)
    }
    attr(psi.mat, "free.pars") <- free.pars
    attr(psi.mat, "free.pars.pos") <- free.pars.pos
    psi.mat
}

## Coefficients between latent endogenous and observed
LAMBDA.Y <- function(model) {
    vars <- attr(model, "vars")
    var.names <- rownames(vars)
    observed <- !vars[,"lat"] 
    latents <- vars[,"lat"] | (!vars[,"lat"] & (vars[,"exo"] |
                                                vars[,"endo"]))
    loadings.y <- lapply(model, `[`, model$op == "=~" & !model$rhs %in%
                         var.names[vars[,"lat"]])
    lambda.y.mat <- matrix(0, ncol = sum(latents), nrow =
                           sum(observed))
    colnames(lambda.y.mat) <- var.names[latents]
    rownames(lambda.y.mat) <- var.names[observed]
    free.pars.pos <- NULL
    free.pars <- loadings.y$free[loadings.y$free > 0L]
    index.matrix <- lambda.y.mat
    index.matrix[] <- seq_along(index.matrix)
    latenterised <- var.names[!vars[,"lat"] & (vars[,"exo"] |
                                                vars[,"endo"])]
    for (i in latenterised) {
        lambda.y.mat[i, i] <- 1
    }
    for (i in seq_len(length(loadings.y$lhs))) {
        lhs <- loadings.y$lhs[i]
        rhs <- loadings.y$rhs[i]
        start.value <- loadings.y$ustart[i]
        is.free <- loadings.y$free[i] > 0L
        lambda.y.mat[rhs, lhs] <-
            loadings.y$ustart[i]
        if (is.free) {
            free.pars.pos <- c(free.pars.pos, index.matrix[rhs, lhs])
        }
    }
    if (length(lambda.y.mat) == 0L) {
        return(II(model))
        free.pars.pos <- NULL
        free.pars <- NULL
    }
    attr(lambda.y.mat, "free.pars") <- free.pars
    attr(lambda.y.mat, "free.pars.pos") <- free.pars.pos
    lambda.y.mat
}

## Variance/covariance of measurement errors of observed endogenous
THETA.EPSILON <- function(model) {
    vars <- attr(model, "vars")
    var.names <- rownames(vars)
    observed <- !vars[,"lat"]
    indicators <- !vars[,"lat"] & (!vars[,"exo"] &
                                                !vars[,"endo"])
    te <- matrix(0, sum(observed), sum(observed))
    obs.names <- var.names[observed]
    ind.names <- var.names[indicators]
    dimnames(te) <- list(obs.names, obs.names)
    var.covar <- lapply(model, `[`, model$op == "~~" & model$lhs
                        %in% ind.names)
    free.pars.pos <- free.pars <- NULL
    index.matrix <- te
    index.matrix[] <- seq_along(index.matrix)
    if (length(te) == 0L) {
        return(matrix(0))
    }
    for (i in seq_len(length(var.covar$lhs))) {
        lhs <- var.covar$lhs[i]
        rhs <- var.covar$rhs[i]
        start.value <- var.covar$ustart[i]
        pos <- var.covar$free[i]
        is.free <- pos > 0L
        if (lhs == rhs) {
            te[lhs, rhs] <- start.value
            if (is.free) {
                free.pars.pos <- c(free.pars.pos, index.matrix[lhs,
                                   rhs])
                free.pars <- c(free.pars, pos)
            }
        } else {
            te[lhs, rhs] <- te[rhs, lhs] <- start.value
            if (is.free) {
                free.pars.pos <- c(free.pars.pos, index.matrix[lhs,
                                   rhs], index.matrix[rhs, lhs])
                free.pars <- c(free.pars, rep(pos, 2))
            }
        }
    }
    attr(te, "free.pars") <- free.pars
    attr(te, "free.pars.pos") <- free.pars.pos
    te
}

## Get the number of degrees of freedom of a given model
model.df <- function(model) {
    ## t-rule
    ## consider equality constraints eventually
    nexo <- length(exogenous(model))
    nvar <- length(observed(model))
    npar <- max(model$free)
    model <- lapply(model, `[`, model$op == "~~" & model$free == 0 & model$lhs == model$rhs)
    nfix <- length(unique(c(model$lhs, model$rhs)))
    ndf <- nvar*(nvar + 1)/2 - npar - nfix*(nfix + 1)/2 
    list(free = npar, 
         fixed = nfix,
         total.df = ndf + npar,
         model.df = ndf)
}

## Function to extract the starting values of the free parameters of a
## model
jumpstart <- function(model) {
    model$ustart[model$free > 0L]
}

## Function to be passed to the optimizer to minimize the log likelihood
iterator <- function(pars, model, sigma, matrep, fun = "fml") {
    model$ustart[model$free > 0L] <- pars
    sigma.hat <- sem.matrices(model, matrep)
    do.call(fun, list(sigma, sigma.hat))
}

## This function takes the 4 basic matrices and returns a covariance
## matrix based on the model parameters.  It can be used to generate a
## 'population' cov matrix for testing purposes or to generate sigmas
## for optimization purposes for instance.

sem.matrices <- function(model, matrep) {
    if (!is.null(model$est)) {
        model$ustart <- model$est
    }
    if (any(is.na(model$ustart))) stop("Please give starting values for
                                       all coefficients in the model.")
    matrep$BE[attr(matrep$BE, "free.pars.pos")] <-
        model$ustart[match(attr(matrep$BE, "free.pars"), model$free)]
    matrep$PS[attr(matrep$PS, "free.pars.pos")] <-
        model$ustart[match(attr(matrep$PS, "free.pars"), model$free)]
    matrep$LY[attr(matrep$LY, "free.pars.pos")] <-
        model$ustart[match(attr(matrep$LY, "free.pars"), model$free)]
    matrep$TE[attr(matrep$TE, "free.pars.pos")] <-
        model$ustart[match(attr(matrep$TE, "free.pars"), model$free)]
    matrep$LY %*% solve(matrep$Inv - matrep$BE) %*% matrep$PS %*%
    solve(t(matrep$Inv - matrep$BE)) %*% t(matrep$LY) + matrep$TE
}

## Maximum-likelihood fitting function.
fml <- function(sigma, sigma.hat) {
    sigma.hat.inv <- try(chol2inv(chol(sigma.hat)), silent = TRUE)
    if (class(sigma.hat.inv) == "try-error") return(Inf)
    res <- (log(det(sigma.hat)) + 
    sum(sigma * sigma.hat.inv) -
    log(det(sigma)) -
    NROW(sigma.hat))
    res
}

## Generalised least squares.
fgls <- function(sigma, sigma.hat) {
    0.5 * sum(diag((solve(sigma) %*% (sigma - sigma.hat))^2))
}

## Function that calls the optimiser and builds a response
#fit.nlm <- function(model, sigma, n, fun = "fml") {
    #model <- lvnfy(model, data = sigma, n = n)
    #pars <- model$ustart[model$free > 0L]
    #SCALE <- ifelse(abs(pars) > 1.00, abs(1/pars), 1)
    #matrep <- matrix.representation(model)
    #rnames <- rownames(sigma.hat(matrep))
    #sigma <- sigma[rnames, rnames]
    #res <- nlm(p = pars, 
                  #f = iterator, 
                  #hessian = TRUE,
                  #check.analyticals = TRUE,
                  #model = model,
                  #fun = fun, 
                  #fscale = SCALE,
                  #sigma = sigma, 
                  #matrep = matrep)
    #hess <- res$hessian
    #std.errors <- se.sem(hess, n)
    #if (!std.errors$trust) warning("Could not compute standard errors.")
    #model$est <- model$ustart
    #model$est[model$free > 0L] <- res$estimate
    #model$se <- rep(0, length(model$id))
    #model$se[model$free > 0L] <- std.errors$values
    #model$z <- model$est / model$se
    #model$z <- ifelse(model$z == Inf, NA, model$z)
    #model$p <- round(do.call("pval", list(model$z)), 3)
    #res$model <- data.frame(model, stringsAsFactors = FALSE)
    #res$chisquared <- n * res$minimum
    #res$df <- model.df(model)$model.df
    #res$n <- n
    #res$pval <- 1 - pchisq(res$chisquared, res$df) 
    #res$estimator <- fun
    #res$sigma.hat <- sem.matrices(model, matrep)
    #res$sigma <- sigma
    #class(res) <- c("sem", "list")
    #invisible(res)
#}

#fit.optim <- function(model, sigma, n, fun = "fml") {
    #model <- lvnfy(model, data = sigma, n = n)
    #pars <- model$ustart[model$free > 0L]
    #SCALE <- ifelse(abs(pars) > 1.00, abs(1.00/pars), 1.00)
    #matrep <- matrix.representation(model)
    #rnames <- rownames(sigma.hat(matrep))
    #sigma <- sigma[rnames, rnames]
    #res <- optim(par = pars, 
                  #fn = iterator, 
                  #gr = gradit,
                  #hessian = TRUE,
                  #model = model,
                  #fun = fun, 
                  #sigma = sigma, 
                  #matrep = matrep,
                  #method = "L-BFGS-B")
    #hess <- res$hessian
    #std.errors <- se.sem(hess, n)
    #model$est <- model$ustart
    #model$est[model$free > 0L] <- res$par
    #model$se <- rep(0, length(model$id))
    #model$se[model$free > 0L] <- std.errors
    #model$z <- model$est / model$se
    #model$z <- ifelse(model$z == Inf, NA, model$z)
    #model$p <- round(do.call("pval", list(model$z)), 3)
    #res$model <- data.frame(model, stringsAsFactors = FALSE)
    #res$chisquared <- n * res$value
    #res$df <- model.df(model)$model.df
    #res$n <- n
    #res$pval <- 1 - pchisq(res$chisquared, res$df) 
    #res$estimator <- fun
    #res$sigma.hat <- sem.matrices(model, matrep)
    #res$sigma <- sigma
    #class(res) <- c("sem", "list")
    #res
#}

#fit.optimx <- function(model, sigma, n, fun = "fml") {
    #model <- lvnfy(model, data = sigma, n = n)
    #pars <- model$ustart[model$free > 0L]
    #SCALE <- ifelse(pars > 1.00, abs(pars/1.00), 1.00)
    #matrep <- matrix.representation(model)
    #rnames <- rownames(sigma.hat(matrep))
    #sigma <- sigma[rnames, rnames]
    #res <- optimx(par = pars, 
                  #fn = iterator, 
                  #gr = gradit,
                  #hessian = TRUE,
                  #model = model,
                  #fun = fun, 
                  #sigma = sigma, 
                  #matrep = matrep,
                  #control = list(factr = 1e-14, all.methods = TRUE),
                  #method = "L-BFGS-B")
    #hess <- res$hessian
    #std.errors <- se.sem(hess, n)
    #model$est <- model$ustart
    #model$est[model$free > 0L] <- res$par
    #model$se <- rep(0, length(model$id))
    #model$se[model$free > 0L] <- std.errors
    #model$z <- model$est / model$se
    #model$z <- ifelse(model$z == Inf, NA, model$z)
    #model$p <- round(do.call("pval", list(model$z)), 3)
    #res$model <- data.frame(model, stringsAsFactors = FALSE)
    #res$chisquared <- n * res$value
    #res$df <- model.df(model)$model.df
    #res$n <- n
    #res$pval <- 1 - pchisq(res$chisquared, res$df) 
    #res$estimator <- fun
    #res$sigma.hat <- sem.matrices(model, matrep)
    #res$sigma <- sigma
    #class(res) <- c("sem", "list")
    #res
#}

## Using stats::nlminb to try and find the best set of parameter
## estimates for a given model. We tried a bunch of other functions, but
## nlminb seems to be the most reliable.
fit.nlminb <- function(model, sigma, n, fun = "fml") {

    ## Setting things up
    model <- lvnfy(model, data = sigma, n = n)
    pars <- model$ustart[model$free > 0L]

    ## Setting the scale (got the idea from lavaan)
    SCALE <- ifelse(abs(pars) > 1.00, abs(1.00/pars), 1)
    matrep <- matrix.representation(model)
    rnames <- rownames(sigma.hat(matrep))

    ## Making sure sigma follows sigma.hat
    sigma <- sigma[rnames, rnames]

    ## Optimizer (finding a set of parameter estimates)
    res <- nlminb(start = pars, 
                  objective = iterator, 
                  gradient = first.derivative.param,
                  ## Using numDeriv::grad works as well
                  #gradient = gradit,
                  model = model,
                  fun = fun, 
                  scale = SCALE,
                  sigma = sigma, 
                  matrep = matrep,
                  control = list(iter.max = 10000, eval.max = 20000))

    ## Calculating the Hessian (standard errors come from this since
    ## this is the inverted of the vcov). There is no reason to expect
    ## this to be pd if there is multiple minima or if the best set of
    ## parameters wasnt found.
    hess <- hessian(func = hessmin, 
                    x = res$par, 
                    sigma = sigma,
                    matrep = matrep,
                    model = model, 
                    fun1 = fun)

    ## Gradients (indicate if the optimizer reach the minima). Here we
    ## used numDeriv::grad, which is a bit slower than
    ## lavaan::first.derivative.param.
    gradd <- grad(func = gradmin, 
                  x = res$par,
                  sigma = sigma,
                  matrep = matrep,
                  model = model,
                  fun = fun)

    ## Checking if any of the eigenvalues are positive.
    ## vcov and standard errors come from the inverse of this matrix. If
    ## it is not invertible, there is no way to calculate the confidence
    ## intervals.
    if (any(eigen(hess)$values < 0L)) {
        stop("Hessian matrix is not positive-definite. \ 
The optimizer did not converge to a unique \
solution. The model is probably empirically \
underidentified.")
    }

    ## Any gradient not close to 0 might indicate that the system did
    ## not reach the minima.
    #if (any(gradd > 1)) { 
        #warning("Some first order derivative are not \
#close to 0. The system probably did \
#not reach the minima. Results might \
#not be thurstworthy.")
    #}

    ## More housecleaning and setting up
    freeid <- model$id[model$free > 0L]
    names(gradd) <- paste(model$lhs[freeid], model$op[freeid],
                         model$rhs[freeid], sep = "")
    colnames(hess) <- rownames(hess) <- names(gradd)

    ## Calculating standard errors
    std.errors <- se.sem(hess, n)
    if (!std.errors$trust) warning("Could not compute standard errors.")

    ## Building the response
    model$est <- model$ustart
    model$est[model$free > 0L] <- res$par
    model$se <- rep(0, length(model$id))
    model$se[model$free > 0L] <- std.errors$values

    ## Calculating z-values
    model$z <- model$est / model$se
    model$z <- ifelse(model$z == Inf, NA, model$z)

    ## Significance of the parameter estimates
    model$p <- round(do.call("pval", list(model$z)), 3)

    out <- model[-1]
    swap <- which(out$op == "~")
    var1 <- out$lhs
    var2 <- out$rhs
    var1[swap] <- out$rhs[swap]
    var2[swap] <- out$lhs[swap]
    out$lhs <- var1
    out$rhs <- var2
    out$op[out$op == "~"] <- "-->"
    out$op[out$op == "~~"] <- "<->"
    out$op[out$op == "=~"] <- "-->"
    out$free <- ifelse(out$free > 0, "free", "fixed")
    out$sig <- rep("*", length(out$lhs))
    out$sig[out$p >= 0.05] <- "ns"
    out$sig[is.na(out$p)] <- NA

    res$model <- data.frame(out, stringsAsFactors = FALSE)
    colnames(res$model) <- c("var1", "op", "var2", "type", "start", "est", "se",
                             "z", "p", "sig")

    ## Calculating MLX2
    res$chisquared <- n * res$objective

    ## Number of degrees of freedom
    res$df <- model.df(model)$model.df

    ## Number of observations
    res$n <- n

    ## Significance of the MLX2
    res$pval <- 1 - pchisq(res$chisquared, res$df) 

    ## Which function was minimised?
    res$estimator <- fun

    ## Calculating the sigma hat with the best set of estimates
    res$sigma.hat <- sem.matrices(model, matrep)

    ## Setting up stuff
    res$sigma <- sigma
    res$grad <- res$gradd
    res$hessian <- hess
    class(res) <- c("sem", "list")
    invisible(res)
}

gradit <- function(pars, model, matrep = matrep, sigma, fun) {
    grad(func = gradmin, x = pars, sigma = sigma,
         model = model, matrep = matrep, fun = fun)
}

gradmin <- function(pars, model, matrep = matrep, sigma, fun) {
    model$ustart[model$free > 0L] <- pars
    sigma.hat <- sem.matrices(model, matrep)
    do.call(fun, list(sigma, sigma.hat))
}

hessmin <- function(pars, model, matrep = matrep, sigma, fun1) {
    model$ustart[model$free > 0L] <- pars
    sigma.hat <- sem.matrices(model, matrep)
    do.call(fun1, list(sigma, sigma.hat))
}

hessit <- function(pars, model, sigma, matrep, fun1) {
    hessian(func = hessmin, x = pars, sigma = sigma,
        model = model, matrep = matrep, fun1 = fun1)
}

se.sem <- function(hess, n) {
    ## From SAS:
    ## http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug
    ## _tcalis_sect072.htm
    hess.inv <- try(chol2inv(chol(hess)))
    if (class(hess.inv) == "try-error") {
        return(list(trust = FALSE, values = rep(NA, NROW(hess))))
    } else {
        return(list(trust = TRUE, values =
                    sqrt(diag(2/(n-1)*hess.inv))))
    }
}

## P-values of parameter estimates
pval <- function(z) {
    2 * (1 - pnorm(abs(z)))
}

## Residuals
residuals.sem <- function(object, ...) {
    round(object$sigma - object$sigma.hat, 4)
}

## Most of this function was taken from lavaan.
RMSEA <- function(fitted.model) {
    rmsea <- sqrt(max(c((fitted.model$chisquared/fitted.model$n)/fitted.model$df - 1/fitted.model$n, 0)))
    lower.lambda <- function(lambda) {
        pchisq(fitted.model$chisquared, df=fitted.model$df, ncp=lambda) - 0.95
    }
    if (fitted.model$df < 1 || lower.lambda(0) < 0) {
        ci.lower <- 0  
    } else {
        lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=fitted.model$chisquared)$root)
        if(inherits(lambda.l, "try-error")) { lambda.l <- NA }
        ci.lower <- sqrt( lambda.l/(fitted.model$n*fitted.model$df) )
    }    
    upper.lambda <- function(lambda) {
        pchisq(fitted.model$chisquared, df=fitted.model$df, ncp=lambda) - 0.05
    }
    N.RMSEA <- max(fitted.model$n, 10000L)
    if (fitted.model$df < 1 || upper.lambda(N.RMSEA) > 0) {
        ci.upper <- 0
    } else {
        lambda.u <- try(uniroot(f=upper.lambda, lower=0, upper=N.RMSEA)$root)
        if(inherits(lambda.u, "try-error")) { lambda.u <- NA }
        ci.upper <- sqrt( lambda.u/(fitted.model$n*fitted.model$df) )
    }
    p <- 1 - pchisq(fitted.model$chisquared, df=fitted.model$df, ncp=(fitted.model$n*fitted.model$df*0.05^2))
    list(RMSEA = rmsea, ci = c(lower = ci.lower, upper = ci.upper), p = p)
}

## Borrowed from lavaan
first.derivative.param <- function(x, model, matrep, sigma, fun) {
    npar <- length(x)
    h <- 10e-6
    dx <- numeric(npar)

    for(i in 1:npar) {
        x.left <- x.left2 <- x.right <- x.right2 <- x
        x.left[i]  <- x[i] - h; x.left2[i]  <- x[i] - 2*h
        x.right[i] <- x[i] + h; x.right2[i] <- x[i] + 2*h
        fx.left   <- iterator(x.left, model = model, matrep = matrep,
                              sigma = sigma, fun = fun)
        fx.left2  <- iterator(x.left2, model = model, matrep = matrep,
                              sigma = sigma, fun = fun)
        fx.right  <- iterator(x.right, model = model, matrep = matrep,
                              sigma = sigma, fun = fun)
        fx.right2 <- iterator(x.right2, model = model, matrep = matrep,
                              sigma = sigma, fun = fun)
        dx[i] <- (fx.left2 - 8*fx.left + 8*fx.right - fx.right2)/(12*h)
    }
    dx
}

summary.sem <- function(object, ...) {
    space <- "\n\n"
    if (object$convergence != 0) {
        txt.con <- "The optimizer did not converge after %s iterations. Results are not trustworthy."
        cat(space, sprintf(txt.con, object$iterations), space)
    } else {
        txt.con <- "The optimizer successfully converged after %s iterations"

        cat(sprintf(txt.con, object$iterations), space)

        cat("MODEL FIT", space)

        cat(sprintf("Method %30s", toupper(object$estimator)), "\n")
        cat(sprintf("Model chi-square %20.3f", object$chisquared), "\n")
        cat(sprintf("P-value %29.3f", object$pval), "\n")
        cat(sprintf("Degrees of freedom %18i", object$df), "\n")
        cat(sprintf("Number of observations %14i", object$n), space)

        cat("PARAMETER ESTIMATES\n\n")

        #cat(sprintf("%30s %10s %10s %10s", "est", "std.error", "z",
                    #"p-value"), "\n")
        #for (i in object$model$id) {
            #cat(sprintf("%-5s %2s %5s %15.3f %10.3f %10.3f %10.3f", object$model$lhs[i],
                        #object$model$op[i], object$model$rhs[i], object$model$est[i],
                        #object$model$se[i], object$model$z[i], object$model$p[i]),
                #"\n")
            #ops <- unique(object$model$op[i:(i+1)])
            #if (length(ops) > 1) {
                #cat("\n")
            #}
        #}
        print(object$model, digits = 3)
    }
}

variables <- function(model) {
    allvars <- unique(c(model$lhs, model$rhs))
    latent <- unique(model$lhs[model$op == "=~"])
    indicators <- unique(c(model$rhs[model$op == "=~" & !model$rhs %in%
                latent]))
    regressions <- lapply(model, `[`, model$op == "~" | (model$op == "=~"
                & model$rhs %in%
                latent))
    exogenous <- unique(regressions$rhs[!regressions$rhs %in%
            regressions$lhs])
    endogenous <- unique(regressions$lhs)
    latent.allvars <- allvars %in% latent
    indicators.allvars <- allvars %in% indicators
    exogenous.allvars <- allvars %in% exogenous
    endogenous.allvars <- allvars %in% endogenous
    #if (!all(exogenous.allvars) & ) exogenous.allvars <- latent.allvars
    #latent.allvars
    out <- cbind(latent.allvars, indicators.allvars, exogenous.allvars,
        endogenous.allvars)
    rownames(out) <- allvars
    colnames(out) <- c("lat", "ind", "exo", "endo")
    out
}

matrix.representation <- function(model) {
    BE <- BETA(model)
    PS <- PSI(model)
    LY <- LAMBDA.Y(model)
    TE <- THETA.EPSILON(model)
    Inv <- II(model)
    list(BE = BE, PS = PS, LY = LY, TE = TE, Inv = Inv)
}

sigma.hat <- function(matrep) {
    matrep$LY %*% solve(matrep$Inv - matrep$BE) %*% matrep$PS %*%
    solve(t(matrep$Inv - matrep$BE)) %*% t(matrep$LY) + matrep$TE
}

boot.lavaan <- function(fitted.model, n) {
    sigma.hat <- fitted.values(fitted.model)$cov
    nobs <- fitted.model@Sample@nobs[[1]]
    mat <- mvrnorm(Sigma = sigma.hat, n = nobs, mu = rep(0,
        NCOL(sigma.hat)), empirical = TRUE)
    est.coefs <- coef(fitted.model)
    X2 <- numeric(n)
    coefs <- matrix(0, nrow = length(est.coefs), ncol = n)
    for (i in 1:n) {
        success <- FALSE
        while (!success) {
            Z <- mat[sample(1:NROW(mat), replace = TRUE), ]
            new.fit <- update(fitted.model, sample.cov =
                                              cov(Z), warn = FALSE)
            if (sum(new.fit@Fit@se > 0) == length(est.coefs)) {
                success <- TRUE
                X2[i] <- new.fit@Fit@test[[1]]$stat
                coefs[, i] <- coef(new.fit)
            }
        }
    }
    p.X2 <- sum(X2 > fitted.model@Fit@test[[1]]$stat)/n
    se <- apply(coefs, 1, sd)
    z <- est.coefs/se
    p.est <- pval(z)
    invisible(list(X2 = X2, p.X2 = p.X2, est = est.coefs, se = se, z = z, p.est =
                   p.est, coefs = coefs))
}
