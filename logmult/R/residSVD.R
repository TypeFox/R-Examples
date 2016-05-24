# Functions to get starting values from residuals obtained from a base model
# They improve residSVD() from gnm:
# 1) By allowing eigen value decomposition in addition to singular value decomposition,
# and support skew-symmetric decomposition
# 2) By supporting different specification of layer effects
# 3) By using a formula more adequate for association models: instead of using working residuals
# that amount to (observed - fitted)/fitted, which are most suited to correlation models,
# decompose log(observed/fitted) by default ('what' argument)
# 4) By filling missing values using an
# Expectation-Maximization (EM) algorithm, as described e.g. by:
# Srebro & Jaakkola, Weighted Low-Rank Approximations, Proceedings of the
# Twentieth International Conference on Machine Learning (ICML-2003),
# Washington DC, 2003.

getResid <- function(model, layer=NULL, what=c("ratio", "residuals"),
                     symmetry=c("asymmetric", "symmetric", "skew-symmetric")) {
  what <- match.arg(what)
  symmetry <- match.arg(symmetry)

  if(what == "residuals") {
      res <- residuals(model, "working")

      weights <- res
      # weights() returns NA for NA cells, which (contrary to model$weights)
      # ensures it's the same lenth as obs
      weights[] <- if(!is.null(model$weights)) as.vector(weights(model, "working")) else 1
      res <- res * weights

      if(!is.null(layer)) {
          res <- res[,, layer]
          weights <- weights[,, layer]
      }
      else if(is.array(res) && dim(res) > 2) {
          res <- apply(res, 1:2, sum, na.rm=TRUE)
          weights <- apply(weights, 1:2, sum, na.rm=TRUE)
      }
      else {
          res <- as.matrix(res)
          weights <- as.matrix(weights)
      }

      res <- res/weights

      if(symmetry == "symmetric")
          res <- (res + t(res))/2
      else if(symmetry == "skew-symmetric")
          res <- res - (res + t(res))/2
  }
  else {
      obs <- model$data
      fitted <- fitted(model)

      weights <- obs
      # weights() returns NA for NA cells, which (contrary to model$weights)
      # ensures it's the same lenth as obs
      weights[] <- if(!is.null(model$weights)) as.vector(weights(model, "working")) else 1
      obs <- obs * weights
      fitted <- fitted * weights

      if(!is.null(layer)) {
          obs <- obs[,, layer]
          fitted <- fitted[,, layer]
          weights <- weights[,, layer]
      }
      else if(is.array(obs) && dim(obs) > 2) {
          obs <- apply(obs, 1:2, sum, na.rm=TRUE)
          fitted <- apply(fitted, 1:2, sum, na.rm=TRUE)
          weights <- apply(weights, 1:2, sum, na.rm=TRUE)
      }
      else {
          obs <- as.matrix(obs)
          fitted <- as.matrix(fitted)
          weights <- as.matrix(weights)
      }

      # This is better than putting an arbitrary large (absolute) value,
      # since in the model the fitted value will not generally be 0,
      # thus will most probably be 1 or more: this is the less disturbing value
      obs[obs == 0] <- 1
      suppressWarnings(res <- log(obs/fitted))

      # Weight SVD/EVD by expected frequencies to prevent cells at the intersection
      # of rows and columns with small counts from influencing too much the result.
      rp <- apply(obs, 1, sum, na.rm=TRUE)/sum(obs, na.rm=TRUE)
      cp <- apply(obs, 2, sum, na.rm=TRUE)/sum(obs, na.rm=TRUE)
      res <- res * sqrt(rp %o% cp)

      if(symmetry == "symmetric")
          res <- (res * weights + t(res * weights))/(weights + t(weights))
      else if(symmetry == "skew-symmetric")
          res <- (res * weights - (res * weights + t(res * weights))/2)/(weights + t(weights))
  }

  res
}

residSVD <- function(model, d=1, layer=NULL, what=c("ratio", "residuals"),
                     tolerance=1e-6, iterMax=500) {
    what <- match.arg(what)

    res <- getResid(model, layer, what)

    miss <- !is.finite(res) | abs(res) < tolerance
    # Make rowMeans() work
    res[is.infinite(res)] <- NA
    res[miss] <- (rowMeans(res, na.rm=TRUE)[row(res)[miss]] + colMeans(res, na.rm=TRUE)[col(res)[miss]])/2

    # Last resort, happens when a full row/column is NA (makes sense with rcL)
    res[!is.finite(res)] <- 0

    if(any(miss)) {
        # Arbitrary value so that convergence test works the first time
        recons.old <- 1

        for(i in 1:iterMax) {
            sv <- svd(res, max(d, nrow(res) - i), max(d, ncol(res) - i))

            # Trick suggested on p. 5 : start with high rank, and decrease it until reaching d
            k <- max(d, length(sv$d) - i)

            recons <- sv$u[, 1:k] %*% (sv$d[1:k] * t(sv$v[, 1:k]))

            res[miss] <- recons[miss]

            if(all(abs(recons.old - recons) < tolerance))
                break

            recons.old <- recons
        }
    }

    sv <- svd(res, d, d)

    result <- t(cbind(sqrt(sv$d[1:d]) * t(sv$u[, 1:d]), sqrt(sv$d[1:d]) * t(sv$v[, 1:d])))
    rownames(result) <- c(rownames(res), colnames(res))
    colnames(result) <- 1:d
    result
}

residEVD <- function(model, d=1, layer=NULL, what=c("ratio", "residuals"),
                     skew=FALSE, tolerance=1e-6, iterMax=500) {
    what <- match.arg(what)

    res <- getResid(model, layer, what, if(skew) "skew" else "symm")
    stopifnot(nrow(res) == ncol(res))

    miss <- !is.finite(res) | abs(res) < tolerance
    # Make rowMeans() work
    res[is.infinite(res)] <- NA
    res[miss] <- rowMeans(res, na.rm=TRUE)[row(res)[miss]]

    # Last resort, happens when a full row/column is NA (makes sense with rcL)
    res[!is.finite(res)] <- 0

    # EVD of skew-symmetric matrices, cf. Constantine & Gower,
    # Graphical Representation of Asymmetric Matrices,
    # J. of the Royal Stat. Soc. Series C (Applied Statistics), 1978
    if(skew) {
        d <- 2 * d

        res <- tcrossprod(res) # This is NN'
    }

    if(any(miss)) {
        # Arbitrary value so that convergence test works the first time
        recons.old <- 1

        for(i in 1:iterMax) {
            eig <- eigen(res, symmetric=TRUE)

            # Trick suggested on p. 5: start with high rank, and decrease it until reaching d
            k <- max(d, length(eig$values) - i)

            recons <- eig$vectors[, 1:k] %*% (eig$values[1:k] * t(eig$vectors[, 1:k]))

            res[miss] <- recons[miss]

            if(all(abs(recons.old - recons) < tolerance))
                break

            recons.old <- recons
        }
    }

    eig <- eigen(res, symmetric=TRUE)

    if(skew)
        result <- t(sign(eig$values[1:d])^(2:1) * sqrt(abs(eig$values[1:d])) * t(eig$vectors[, 1:d]))
    else
        result <- t(sqrt(abs(eig$values[1:d])) * t(eig$vectors[, 1:d]))
    rownames(result) <- rownames(res)
    colnames(result) <- 1:d
    result
}

residSVDL <- function(model, d=1, layer.effect=c("homogeneous.scores", "heterogeneous", "none"),
                      tolerance=1e-6, iterMax=500) {
    layer.effect <- match.arg(layer.effect)
    res <- residuals(model, "working")

    stopifnot(is.array(res) & length(dim(res)) == 3)

    if(layer.effect == "heterogeneous") {
        result <- apply(simplify2array(lapply(1:dim(res)[3], function(i)
                      residSVD(model, d=d, layer=i, tolerance=tolerance, iterMax=iterMax))), 2, cbind)
        rownames(result) <- outer(c(paste(names(dimnames(res))[1], rownames(res), sep=""),
                                    paste(names(dimnames(res))[2], colnames(res), sep="")),
                                  paste(names(dimnames(res))[3], dimnames(res)[[3]], sep=""),
                                  paste, sep=":")

        # In the model we group scores by row/column, not by layer
        result <- result[(seq(dim(res)[3]) - 1) * (nrow(res) + ncol(res)) +
                          rep(seq(nrow(res) + ncol(res)), each=dim(res)[3]),]
    }
    else {
        result <- residSVD(model, d)
    }

    result
}

residEVDL <- function(model, d=1, layer.effect=c("homogeneous.scores", "heterogeneous", "none"),
                      skew=FALSE, tolerance=1e-6, iterMax=500) {
    layer.effect <- match.arg(layer.effect)
    res <- residuals(model, "working")

    stopifnot(is.array(res) & length(dim(res)) == 3)

    if(layer.effect == "heterogeneous") {
        result <- apply(simplify2array(lapply(1:dim(res)[3], function(i)
                      residEVD(model, d=d, layer=i, skew=skew, tolerance=tolerance, iterMax=iterMax))), 2, cbind)
        rownames(result) <- outer(paste(names(dimnames(res))[1], "|", names(dimnames(res))[2],
                                        rownames(res), sep=""),
                                  paste(names(dimnames(res))[3], dimnames(res)[[3]], sep=""),
                                  paste, sep=":")

        # In the model we group scores by row/column, not by layer
        result <- result[(seq(dim(res)[3]) - 1) * nrow(res) +
                          rep(seq(nrow(res)), each=dim(res)[3]),]
    }
    else {
        result <- residEVD(model, d, skew=skew)
    }

    result
}
