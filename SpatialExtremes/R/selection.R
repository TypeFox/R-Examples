TIC.default <- function(object, ..., k = 2){

  all.objects <- c(list(object), list(...))
  n.models <- length(all.objects)
  tic <- rep(NA, n.models)

  for (i in 1:n.models){

    object <- all.objects[[i]]
    ihessian <- object$ihessian
    var.score <- object$var.score

    if (!is.null(ihessian) && !is.null(var.score)){
      penalty <- var.score %*% ihessian
      tic[i] <- deviance(object) + k * sum(diag(penalty))
    }

  }

  names(tic) <- as.character(sys.call())[2:(n.models + 1)]
  tic <- tic[order(tic)]

  return(tic)
}

anova.maxstab <- function(object, object2, method = "RJ",
                          square = "chol", ...){

  if (is.null(object$var.cov) || is.null(object2$var.cov))
    return("Standard errors must be available for both fitted models.")

  if (!(method %in% c("RJ", "CB")))
    stop("'method' must be one of 'RJ' or 'CB'")

  if (!(square %in% c("chol", "svd")))
    stop("'square' must be one of 'chol' or 'svd'")

  ##Check if object and object2 are nested
  n.estim <- length(fitted(object))
  n.estim2 <- length(fitted(object2))

  if (object$cov.mod != object2$cov.mod)
    stop("Models are not nested.")

  if (n.estim == n.estim2)
    stop("Models are not nested.")

  else{
    if (n.estim > n.estim2){
      M0 <- object2
      M1 <- object
      model0 <- deparse(substitute(object2))
      model1 <- deparse(substitute(object))
    }

    else{
      M0 <- object
      M1 <- object2
      model0 <- deparse(substitute(object))
      model1 <- deparse(substitute(object2))
    }
  }

  models <- c(model0, model1)

  ihessian <- M1$ihessian
  var.cov <- M1$var.cov

  ## Beware a isotropic Smith model has a param named "cov" and not
  ## "covXX", need to fix it in such cases
  if (M0$iso != M1$iso){
    if (M0$iso){
      idx <- which(names(M0$fitted.values) == "cov")
      names(M0$fitted.values)[idx] <- "cov11"

      M0$fixed <- c(M0$fixed, cov12 = 0, M0$par["cov22"])
    }

    if (M1$iso){
      idx <- which(names(M1$fitted.values) == "cov")
      names(M1$fitted.values)[idx] <- "cov11"

      M1$fixed <- c(M1$fixed, cov12 = 0, M1$par["cov22"])
    }
  }
  removed.param <- which(!(names(M1$fitted.values) %in%
                             names(M0$fitted.values)))
  removed.param <- names(M1$fitted.values)[removed.param]

  if (method == "RJ")
    Dev <- c(deviance(M0), deviance(M1))

  else{
    theta0 <- M0$param
    theta0 <- theta0[colnames(ihessian)]

    jac <- M1$var.score
    ijac <- try(solve(jac), silent = TRUE)
    hessian <- try(solve(ihessian), silent = TRUE)

    if (is.matrix(ijac) && is.matrix(hessian)){

      ivar.cov <- hessian %*% ijac %*% hessian

      if (square == "svd"){
        svd.hessian <- svd(hessian)
        svd.ivar.cov <- svd(ivar.cov)
        M <- svd.hessian$u %*% diag(sqrt(svd.hessian$d)) %*%
          t(svd.hessian$u)
        Madj <- svd.ivar.cov$u %*% diag(sqrt(svd.ivar.cov$d)) %*%
          t(svd.ivar.cov$u)
      }

      else{
        M <- try(chol(hessian), silent = TRUE)
        Madj <- try(chol(ivar.cov), silent = TRUE)
      }

      if (is.matrix(M) && is.matrix(Madj)){
        C <- solve(M) %*% Madj
        colnames(C) <- rownames(C) <- colnames(ihessian)

        theta0.adj <- as.numeric(fitted(M1) + C %*% (theta0 - fitted(M1)))
        names(theta0.adj) <- colnames(C)
        theta0.adj <- c(list(p = theta0.adj), M1$fixed)
        Dev <- c(2 * do.call(M1$nllh, theta0.adj), deviance(M1))

        c <- (M1$param[removed.param] - M0$fixed[removed.param])  %*%
          solve(var.cov[removed.param, removed.param]) %*%
            (M1$param[removed.param] - M0$fixed[removed.param]) /
              ((fitted(M1) - theta0) %*% ivar.cov %*%
             (fitted(M1) - theta0))
        Dev <- Dev * c

      }

      else{
        warning("Impossible to get the Cholesky or the svd decomposition")
        Dev <- c(NA, deviance(M1))
      }
    }

    else{
      warning("Matrices H or/and J are singular")
      Dev <- c(NA, deviance(M1))
    }
  }

  diffDev <- Dev[1] - Dev[2]
  MDf <- c(length(fitted(M0)), length(fitted(M1)))
  Df <- diff(MDf)

  if (method == "RJ"){
    ##Under model misspecification, the distribubion of
    ##diffDev is \sum_i \eigen_i chi^2_1
    ihessian <- M1$ihessian
    hessian <- solve(ihessian[removed.param,removed.param])
    var.cov <- var.cov[removed.param,removed.param]

    Q <- var.cov %*% hessian
    eigen.val <- eigen(Q)$values

    pvalue <- pchisq(Df * diffDev / sum(eigen.val), df = Df, lower.tail = FALSE)
  }

  else{
    eigen.val <- rep(1, Df)
    pvalue <- pchisq(diffDev, Df, lower.tail = FALSE)
  }


  table <- data.frame(MDf, Dev, c(NA, Df), c(NA, diffDev),
                      c(NA, pvalue))

  dimnames(table) <- list(models, c("MDf", "Deviance", "Df",
                                    "Chisq", "Pr(> sum lambda Chisq)"))

  if (method == "RJ")
    structure(table, heading = paste(c("Eigenvalue(s):", round(eigen.val, 2),
                       "\n\nAnalysis of Variance Table"), collapse = " "),
              class = c("anova", "data.frame"))

  else
    structure(table, heading = "Analysis of Variance Table",
              class = c("anova", "data.frame"))
}


anova.spatgev <- function(object, object2, method = "RJ",
                          square = "chol", ...){

  if (is.null(object$var.cov) || is.null(object2$var.cov))
    return("Standard errors must be available for both fitted models.")

  if (!(method %in% c("RJ", "CB")))
    stop("'method' must be one of 'RJ' or 'CB'")

  if (!(square %in% c("chol", "svd")))
    stop("'square' must be one of 'chol' or 'svd'")

  ##Check if object and object2 are nested
  n.estim <- length(fitted(object))
  n.estim2 <- length(fitted(object2))

  if (n.estim == n.estim2)
    stop("Models are not nested.")

  else{
    if (n.estim > n.estim2){
      M0 <- object2
      M1 <- object
      model0 <- deparse(substitute(object2))
      model1 <- deparse(substitute(object))
    }

    else{
      M0 <- object
      M1 <- object2
      model0 <- deparse(substitute(object))
      model1 <- deparse(substitute(object2))
    }
  }

  models <- c(model0, model1)

  ihessian <- M1$ihessian
  var.cov <- M1$var.cov

  removed.param <- which(!(names(M1$fitted.values) %in%
                             names(M0$fitted.values)))
  removed.param <- names(M1$fitted.values)[removed.param]

  if (method == "RJ")
    Dev <- c(deviance(M0), deviance(M1))

  else{
    theta0 <- M0$param
    theta0 <- theta0[colnames(ihessian)]

    jac <- M1$var.score
    ijac <- try(solve(jac), silent = TRUE)
    hessian <- try(solve(ihessian), silent = TRUE)

    if (is.matrix(ijac) && is.matrix(hessian)){

      ivar.cov <- hessian %*% ijac %*% hessian

      if (square == "svd"){
        svd.hessian <- svd(hessian)
        svd.ivar.cov <- svd(ivar.cov)
        M <- svd.hessian$u %*% diag(sqrt(svd.hessian$d)) %*%
          t(svd.hessian$u)
        Madj <- svd.ivar.cov$u %*% diag(sqrt(svd.ivar.cov$d)) %*%
          t(svd.ivar.cov$u)
      }

      else{
        M <- try(chol(hessian), silent = TRUE)
        Madj <- try(chol(ivar.cov), silent = TRUE)
      }

      if (is.matrix(M) && is.matrix(Madj)){
        C <- solve(M) %*% Madj
        colnames(C) <- rownames(C) <- colnames(ihessian)

        theta0.adj <- as.numeric(fitted(M1) + C %*% (theta0 - fitted(M1)))
        names(theta0.adj) <- colnames(C)
        theta0.adj <- c(list(p = theta0.adj), M1$fixed)
        Dev <- c(2 * do.call(M1$nllh, theta0.adj), deviance(M1))

        c <- (M1$param[removed.param] - M0$fixed[removed.param])  %*%
          solve(var.cov[removed.param, removed.param]) %*%
            (M1$param[removed.param] - M0$fixed[removed.param]) /
              ((fitted(M1) - theta0) %*% ivar.cov %*%
             (fitted(M1) - theta0))
        Dev <- Dev * c

      }

      else{
        warning("Impossible to get the Cholesky or the svd decomposition")
        Dev <- c(NA, deviance(M1))
      }
    }

    else{
      warning("Matrices H or/and J are singular")
      Dev <- c(NA, deviance(M1))
    }
  }

  diffDev <- Dev[1] - Dev[2]
  MDf <- c(length(fitted(M0)), length(fitted(M1)))
  Df <- diff(MDf)

  if (method == "RJ"){
    ##Under model misspecification, the distribubion of
    ##diffDev is \sum_i \eigen_i chi^2_1
    ihessian <- M1$ihessian
    hessian <- solve(ihessian[removed.param,removed.param])
    var.cov <- var.cov[removed.param,removed.param]

    Q <- var.cov %*% hessian
    eigen.val <- eigen(Q)$values

    pvalue <- pchisq(Df * diffDev / sum(eigen.val), df = Df, lower.tail = FALSE)
  }

  else{
    eigen.val <- rep(1, Df)
    pvalue <- pchisq(diffDev, Df, lower.tail = FALSE)
  }


  table <- data.frame(MDf, Dev, c(NA, Df), c(NA, diffDev),
                      c(NA, pvalue))

  dimnames(table) <- list(models, c("MDf", "Deviance", "Df",
                                    "Chisq", "Pr(> sum lambda Chisq)"))

  if (method == "RJ")
    structure(table, heading = paste(c("Eigenvalue(s):", round(eigen.val, 2),
                         "\n\nAnalysis of Variance Table"), collapse = " "),
              class = c("anova", "data.frame"))

  else
    structure(table, heading = "Analysis of Variance Table",
              class = c("anova", "data.frame"))
}

