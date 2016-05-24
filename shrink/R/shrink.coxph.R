#' Global, Parameterwise and Joint Shrinkage of Regression Coefficients for Fit Object's of Class \code{coxph}.
#'
#' Internal function for \code{shrink}.
#'
#' @inheritParams shrink
#'
#' @return A list with the following elements
#' \tabular{ll}{
#'     ShrinkageFactors \tab cdot \cr
#'     ShrinkageFactorsVCOV \tab sfit$var \cr
#'     ShrunkenRegCoef \tab cdot * beta \cr
#'     postfit \tab pfit \cr
#' }
#'
#' @keywords internal
#' @export
shrink.coxph <-
  function(fit, type, method, join, notes, postfit)
{
  beta <- coef(fit)
  n <- fit$n
  if (is.null(mwts <- weights(fit))) mwts <- rep(1L, n)
  mmat <- model.matrix(fit)
  y <- fit$y
  mmeth <- fit$method
  mcont <- coxph.control()
  #  if (inherits(fit, "mfp")) varnames <- names(fit$fit$coefficients) else
  varnames <- names(fit$coefficients)

  if (!is.null(join)) {                                                           # for rcs
    join2 <- sapply(join, function(x) sub(x = x, pattern = "rcs(", replacement = "",   # extract names of rcs variables
                                          fixed = TRUE), simplify = FALSE)
    for (i in 1L:length(join)) {
      comp <- join2[[i]] == join[[i]]
      join2[[i]][!comp] <- sub(join2[[i]][!comp], pattern = ")", replacement = "", fixed = TRUE)

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
        #          if (pre != strsplit(strsplit(varnames[index[i]], ")")[[1L]][2L], "'")[[1L]][1L]) {
        #            pre <- strsplit(varnames[index[i]], ")")[[1L]][2L]
        #          }
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
        stop("variables in 'fit' and variables listed in 'join' do not match. Note: For objects of class 'mfp' use variable names from 'names(coef(fit))'.")

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
        stop("variables in 'fit' and variables listed in 'join' do not match. Note: For objects of class 'mfp' use variable names from 'names(coef(fit))'.")
    }
  }

  if (method == "dfbeta") {
    beta.hat <- matrix(rep(beta, each = n), nrow = n) - cbind(residuals(fit, "dfbeta"), deparse.level = 0L)
    # works without matrix, but then the results are slightly different
  } else {
    beta.hat <- array(dim = c(n, length(beta)), dimnames = list(1L:n, names(beta)))
    for (i in 1L:n) beta.hat[i, ] <-  coxph.fit(y = y[-i], x = mmat[-i, , drop = FALSE],
                                                strata = NULL, control = mcont,
                                                method = mmeth, rownames = NULL)$coefficients
    #    coxph(fit$y[-i,] ~ fit$x[-i,])$coefficients
  }

  eta <- if (type == "global") { cbind(rowSums(mmat * beta.hat), deparse.level = 0L) } else {
    if (is.null(join))  { mmat * beta.hat } else {
      shrinkage <- matrix(NA, nrow = n, ncol = length(varnames) - length(unlist(join)) + length(join))
      Joined <- rep(0L, length(varnames))
      index <- rep(0L, times = length(varnames))
      names(index) <- varnames
      prod <- (mmat * beta.hat)

      for (i in 1L:length(join)) {
        joined <- varnames %in% join[[i]]
        index[joined] <- i
        Joined <- Joined + joined
        shrinkage[,i] <- rowSums(prod[, joined, drop = FALSE])
      }
      if (all(varnames %in% unlist(join)))
        dimnames(shrinkage) <- list(1L:nrow(shrinkage), paste("join", unlist(lapply(join, function(X) X[1L])), sep="."))
      else {                        # shrinkage factors for variables not listed in join
        shrinkage[, -c(1L:length(join))] <- prod[, !Joined, drop = FALSE]
        dimnames(shrinkage)[[2L]] <- c(paste("join", unlist(lapply(join, function(X) X[1L])), sep="."), varnames[!Joined])
      }
      shrinkage
    }
  }

  cdot <- (sfit <- coxph.fit(y = y, x = eta, strata = NULL, control = mcont, method = mmeth,
                             rownames = NULL))$coefficients
  #  cdot <- (sfit <- coxph(y ~ ., data = data.frame(eta)))$coefficient
  if (type == "global") dimnames(sfit$var) <- list("global", "global") else dimnames(sfit$var) <- list(dimnames(eta)[[2]], dimnames(eta)[[2]])

  if (type == "parameterwise" && !is.null(join)) {                                # sort res
    cdotsort2 <- c(varnames[!Joined], paste("join", rep(unlist(lapply(join, function(x) x[1L])),
                                                        unlist(lapply(join, function(x) length(x)))), sep="."))
    cdotsort2 <- cdot[match(cdotsort2, names(cdot))]
    #    if ((sum(Joined == 0L)) == 0L) { names(cdotsort2) <- unlist(join) } else {
    if (all(Joined)) names(cdotsort2) <- unlist(join) else
      names(cdotsort2)[-c(0L:sum(Joined == 0L))] <- unlist(join)
      cdot <- cdotsort2[varnames]
  }

  if (postfit) {
    pfit <- if (type == "global" || !is.null(join))
      coxph(y ~ ., data = data.frame(mmat / rep(cdot, each = n))) else
      coxph(y ~ ., data = data.frame(mmat * beta.hat / rep(beta, each = n)))
  } else pfit <- NULL

  res <- list(ShrinkageFactors = cdot, ShrinkageFactorsVCOV = sfit$var,
              ShrunkenRegCoef = cdot * beta, postfit = pfit)
  return(res)
}
