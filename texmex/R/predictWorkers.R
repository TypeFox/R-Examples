addCov <- function(res, X){ # used in linearPredictors.* to add covariates to columns reported in output
  if(!is.null(dim(X))){
    if(dim(X)[2] > 1){
       cov <- X[,colnames(X) != "(Intercept)", drop=FALSE]
       res <- cbind(res, cov)
       if(is.vector(cov)) colnames(res)[dim(res)[2]] <- colnames(X)[colnames(X) != "(Intercept)"]
    }
    else {
      if( any(X != 1) ){
        res <- cbind(res, X)
      }
    }
  }
  res
}

texmexMakeParams <-
    # Take parameter vector and list of datasets to compute
    # parameters for each row of the data.
function(co, data){

    np <- length(data)
    p <- vector('list', length=np)

    wh <- 1
    for (i in 1:np){
        which <- wh:(wh -1 + ncol(data[[i]]))
        p[[i]] <- c(co[which] %*% t(data[[i]]))
        wh <- wh + ncol(data[[i]])
    }

    do.call('cbind', p)
}

texmexMakeCovariance <-
    # Get covariance matrix for each parameter and get coefficents for each row of the data
function(cov, data){
    covs <- v <- vector('list', length=length(data))

    # First get covariance matrix for each main parameter (e.g. phi=a'X1, xi=b'X2)
    wh <- 1
    for (i in 1:length(covs)){
        which <- wh:(wh -1 + ncol(data[[i]]))
        covs[[i]] <- as.matrix(cov[which, which])
        # covs[[]i] contains the block of the full covariance which relates to parameter[i]
 
        # Get the variance of the linear predictors
        v[[i]] <- rowSums((data[[i]] %*% covs[[i]]) * data[[i]])

        wh <- wh + ncol(data[[i]])
    }
    names(v) <- names(D)
    # Each element of v contains the variance for one linear predictor for every observation

    # We now need the off-diagonal elements of the covariance matrix. The dimensions
    # of the covariance will depend on the length of object$data$D

    getOffDiagonal <- function(k, x1, x2, cov){
        covar <- 0
        for (i in 1:ncol(x1)){
            for (j in 1:ncol(x2)){
                covar <- covar + x1[k, i] * x2[k, j] * cov[i, ncol(x1) + j]
            } # Close for j
        } # Close for i
        covar
    } # Close getOffDiagonal

    getCovEntry <- function(data){
        # Recursively produce off-diagonal elements of the covariance.
        # These are computed by row.
        x1 <- data[[1]]
        data[[1]] <- NULL
        n <- length(data)

        res <- vector('list', length=n)
        for (i in 1:n){
            res[[i]] <- sapply(1:nrow(x1), getOffDiagonal, x1, data[[i]], cov)
        }
        if (length(data) > 1){
            res <- c(res, getCovEntry(data))
        }
        else {
            res
        }
    } # Close getCovEntry

    co <- getCovEntry(data)

    # Now need to restructure to return a list of covariance matrices,
    # one element for every observation. The covariance is of the linear
    # predictors.
    getM <- function(i, variance, covariance){
        va <- lapply(variance, function(x){ x[[i]] })
        co <- lapply(covariance, function(x){ x[[i]] })

        res <- diag(unlist(va))
        res[upper.tri(res)] <-  res[lower.tri(res)] <- unlist(co)
        res
    }

    lapply(1:nrow(data[[1]]), getM, v, co)
}

texmexPredictSE <- function(x){
    # x is an object returned by texmexMakeCovariance. It
    # is a list, each element of which is a covariance matrix.

    se <- sapply(x, function(x){ sqrt(diag(x)) })
}

texmexMakeCI <-
    # Compute CIs from point estimates and standard errors
function(params, ses, alpha){
    z <- qnorm(1 - alpha/2)

    ci <- lapply(1:ncol(params), function(i, x, s, z){
                                     cbind(x[, i] - z*s[, i], x[, i] + z*s[, i])
                                 }, x=params, s=ses, z=z)
    ci <- do.call('cbind', ci)

    nms <- rep('', ncol(ci))
    nms[(1:ncol(ci)) %% 2 == 1] <- paste(colnames(params), '.lo', sep = '')
    nms[(1:ncol(ci)) %% 2 == 0] <- paste(colnames(params), '.hi', sep = '')
    colnames(ci) <- nms

    ci
}

texmexMakeNewdataD <- function(x, newdata){
    if (is.null(newdata)){
        res <- x$data$D
    }
    else {
        xl <- function(i, fo, data, xlev){
                if (length(xlev[[i]]) > 0){
                  model.matrix(fo[[i]], data, xlev=xlev[[i]])
                }
                else {
                  model.matrix(fo[[i]], data)
                }
              } # Close function
        fo <- x$formulae
        res <- lapply(1:length(fo), xl, fo=fo, data=newdata, xlev=x$xlevels)
        names(res) <- names(fo)
    }

    invisible(res)
}

texmexMakeCISim <- function(x, alpha, object, sumfun, M){
    if (is.null(sumfun)){
        sumfun <- function(x){
            c(Mean=mean(x), quantile(x, prob=c(.50, alpha/2,  1 - alpha/2)) )
        }
        neednames <- TRUE
    } else {
        neednames <- FALSE
    }

    t(apply(x, 2, sumfun))
}

