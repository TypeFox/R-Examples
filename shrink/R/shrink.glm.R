#' Global, Parameterwise and Joint Shrinkage of Regression Coefficients for Fit Object's of Class \code{glm}.
#'
#' Internal function for \code{shrink}.
#'
#' @inheritParams shrink
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#'     ShrinkageFactors \tab cdot[rmInt] \cr
#'     ShrinkageFactorsVCOV \tab summary.glm(sfit)$cov.scaled \cr
#'     ShrunkenRegCoef \tab shcoef \cr
#'     postfit \tab pfit \cr
#' }
#'
#' @keywords internal
#' @export
shrink.glm <-
  function(fit, type, method, join, notes, postfit)
{
  if (inherits(fit, "mfp")) {                                           # mfp inconsistency
    if (notes) cat("Note, that in fit objects of class 'mfp' summary(fit) may be != print(fit). shrink uses variable names from 'names(coef(fit))'.\n\n")
    beta <- coef(fit$fit)
    mmat <- model.matrix(fit$fit)
    varnames <- colnames(mmat) <- names(beta) <- names(fit$coefficients)
  } else {
    beta <- coef(fit)
    mmat <- model.matrix(fit)
    varnames <- names(fit$coefficients)
  }

  n <- nobs(fit)
  y <- fit$y
  moff <- fit$offset
  mwts <- weights(fit)
  if (any(mwts==0)) stop("All 'weights' used in fit have to be > 0.")
  mfam <- family(fit)
  #  if (!mfam[[1]]  %in% c("gaussian", "binomial")) stop(paste("family '", mfam[[1L]], "'
  #     is not supported", sep=" "))

  if (!is.null(join)) {                                                           # for rcs
    join2 <- sapply(join, function(x) sub(x = x, pattern = "rcs(", replacement = "",      # extract names of rcs variables
                                          fixed = TRUE), simplify = FALSE)
    for(i in 1L:length(join)) {
      comp <- !(join2[[i]] == join[[i]])
      join2[[i]][comp] <- sub(join2[[i]][comp], pattern = ")", replacement = "", fixed = TRUE)

      if (length(ex <- regexpr(join2[[i]][comp], pattern = ",", fixed = TRUE)) > 0L)       # if rcs with position of knots is specified, if for e.g. "rcsa"
        for(l in 1L:length(ex)) { if (ex[l] > 0L) join2[[i]][comp][l] <- substr(x = join2[[i]][comp][l], start = 1L, stop = ex[l] - 1L) }

      join2[[i]][join[[i]] %in% varnames] <- join[[i]][join[[i]] %in% varnames]   # if too much was removed
    }

    rcs.test <- isrcs(fit)                              # test for wrong non-rcs variables
    if (rcs.test[[1L]] && any(unlist(join) != unlist(join2))) {
      varnames2 <- varnames
      index <- ((1L:length(varnames))[as.vector(rcs.test[[2L]])])
      #      pre <- strsplit(varnames[index[1L]], ")")[[1L]][2L]
      for(i in 1L:length(index))  {
        #        if (i != 1L) {
        #          if (pre != strsplit(strsplit(varnames[index[i]], ")")[[1L]][2L], "'")[[1L]][1L])
        #            pre <- strsplit(varnames[index[i]], ")")[[1L]][2L]
        #        }
        #        varnames2[index[i]] <- paste(strsplit(varnames[index[i]], ")")[[1L]][2L], sep="")
        #        varnames2[index[i]] <- strsplit(varnames2[index[i]], "'")[[1L]][1L]

        varnames2[index[i]] <- sub(x = varnames2[index[i]], pattern = "rcs(", replacement = "", fixed = TRUE)
        if ((ex2 <- regexpr(varnames2[index[i]], pattern = ",", fixed = TRUE)[1L]) > 0L) {
          varnames2[index[i]] <- substr(x = varnames2[index[i]], start = 1L, stop = ex2 - 1L) } else {
            varnames2[index[i]] <- substr(x = varnames2[index[i]], start = 1L,
                                          stop = regexpr(varnames2[index[i]], pattern = ")", fixed = TRUE)[1] - 1L) }
      }

      if (!all(unlist(join2) %in% unique(varnames2)) || !all(unlist(join) %in% attr(terms(fit), "term.labels")))    # use identical rcs call
        stop("variables in 'fit' and variables listed in 'join' do not match. Note: For objects of class 'mfp' use variable names from 'names(coef(fit))'.)")

      rcsch <- lapply(join, function(x) !(x %in% varnames))            # correct rcs names
      if (sum(unlist(rcsch)) >= 1L) {
        ext <- unlist(join)[!(unlist(join) %in% varnames)]
        for (j in length(ext):1L) {
          rcsch2 <- lapply(join, function(x) x == ext[j])
          rcsch3 <- which(unlist(lapply(rcsch2, function(x) sum(x))) == 1L)
          join[[rcsch3]] <- c(join[[rcsch3]][!rcsch2[[rcsch3]]],
                              varnames[grep(pattern = join[[rcsch3]][rcsch2[[rcsch3]]],
                                            x = varnames, fixed = TRUE)])
        }
      }
    } else {
      if (!all(unlist(join) %in% varnames))
        stop("variables in 'fit' and variables listed in 'join' do not match. Note: For objects of class 'mfp' use variable names from 'summary(fit)'.)")
    }
  }

  if (intercept <- attr(terms(fit), "intercept"))  varnames <- varnames[!(varnames %in% c("(Intercept)", "Intercept"))]

  if (method == "dfbeta") {
    if (inherits(fit, "mfp")) {                                         # 1) mfp inconsistency, 2) dfbeta only works with numeric fit$y
      beta.hat <- matrix(rep(beta, each = n), nrow = n) - dfbeta(fit$fit)
    } else {
      beta.hat <- matrix(rep(beta, each = n), nrow = n) - dfbeta(fit) } } else {  # intercept is not considered
        # works without matrix, but then the results are slightly different
        beta.hat <- array(dim = c(n, length(beta)), dimnames = list(1L:n, names(beta)))
        for (i in 1L:n) beta.hat[i, ] <- glm.fit(y = y[-i], x = mmat[-i, , drop = FALSE],
                                                 family = mfam, offset = moff[-i],
                                                 weights = mwts[-i])$coefficients
      }

  beta.hat[, intercept] <- 1L
  rmInt <- if (intercept) -1L else TRUE

  eta <- if (type == "global") {
    j <- c(intercept, rep(intercept + 1L, length(beta) - intercept))
    cbind(if (intercept) 1 else NULL, rowSums((mmat * beta.hat)[, rmInt, drop = FALSE]),
          deparse.level = 0L)
  } else {
    if (is.null(join))  {
      j <- TRUE
      mmat * beta.hat
    } else {
      j <- TRUE
      shrinkage <- matrix(NA, nrow = n, ncol = length(varnames) - length(unlist(join)) + length(join))
      Joined <- rep(0L, length(varnames))
      index <- rep(0L, times = length(varnames))
      names(index) <- varnames
      prod <- (mmat * beta.hat)[, rmInt, drop = FALSE]

      for (i in 1L:length(join)) {
        joined <- varnames %in% join[[i]]
        index[joined] <- i
        Joined <- Joined + joined
        shrinkage[,i] <- rowSums(prod[, joined, drop = FALSE])
      }
      if (all(varnames %in% unlist(join))) {
        dimnames(shrinkage) <- list(1L:nrow(shrinkage), paste("join", unlist(lapply(join, function(X) X[1L])), sep="."))
      } else {                        # shrinkage factors for variables not listed in join
        shrinkage[, -c(1L:length(join))] <- prod[, !Joined, drop = FALSE]
        dimnames(shrinkage)[[2L]] <- c(paste("join", unlist(lapply(join, function(X) X[1L])), sep="."), varnames[!Joined])
      }
      cbind(intercept, shrinkage)
    }
  }

  cdot <- (sfit <- glm.fit(x = eta, y = y, family = mfam))$coefficients
  if (type == "global") names(sfit$coefficients) <- if (intercept) c("intercept", "global") else "global"

  if (type == "parameterwise" && !is.null(join)) {              # sort res, intercept is not shown in the results
    cdotsort <- cdot[rmInt]
    cdotsort2 <- c(varnames[!Joined], paste("join", rep(unlist(lapply(join, function(x) x[1L])),
                                                          unlist(lapply(join, function(x) length(x)))), sep="."))
    cdotsort2 <- cdotsort[match(cdotsort2, names(cdotsort))]
    if (all(Joined)) names(cdotsort2) <- unlist(join) else
      names(cdotsort2)[-c(0L:sum(Joined == 0L))] <- unlist(join)
    if (intercept) cdot <- c(cdot[1L], cdotsort2[varnames]) else cdot <- cdotsort2[varnames]
  }

  shcoef <- if (intercept) {
    beta[1L] <- NA
    shr <- cdot[j] * beta
    lp <- colSums((t(mmat) * shr), na.rm = TRUE)
    shr[1L] <- glm.fit(x = matrix(1L, nrow = n), y = y, offset = lp,
                       family = mfam, weights = mwts)$coefficients
    shr
  } else cdot[j] * beta

  if (!intercept && postfit) {
    pfit <- if (type == "global" || !is.null(join))
      glm(y ~ -1 + ., data = data.frame(mmat / rep(cdot, each = n)), family = mfam) else
        glm(y ~ -1 + ., data = data.frame(mmat * beta.hat / rep(beta, each = n)), family = mfam)
  } else pfit <- NULL

  res <- list(ShrinkageFactors = cdot[rmInt], ShrinkageFactorsVCOV = summary.glm(sfit)$cov.scaled,
              ShrunkenRegCoef = shcoef, postfit = pfit)
  return(res)
}
