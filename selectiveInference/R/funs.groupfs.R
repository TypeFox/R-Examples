#' Select a model with forward stepwise.
#'
#' This function implements forward selection of linear models almost identically to \code{\link[stats]{step}} with \code{direction = "forward"}. The reason this is a separate function from \code{\link{fs}} is that groups of variables (e.g. dummies encoding levels of a categorical variable) must be handled differently in the selective inference framework.
#'
#' @param x Matrix of predictors (n by p).
#' @param y Vector of outcomes (length n).
#' @param index Group membership indicator of length p. Check that \code{sort(unique(index)) = 1:G} where \code{G} is the number of distinct groups.
#' @param maxsteps Maximum number of steps for forward stepwise.
#' @param sigma Estimate of error standard deviation for use in AIC criterion. This determines the relative scale between RSS and the degrees of freedom penalty. Default is NULL corresponding to unknown sigma. When NULL, \code{link{groupfsInf}} performs truncated F inference instead of truncated \eqn{\chi}. See \code{\link[stats]{extractAIC}} for details on the AIC criterion.
#' @param k Multiplier of model size penalty, the default is \code{k = 2} for AIC. Use \code{k = log(n)} for BIC, or \code{k = 2log(p)} for RIC (best for high dimensions, when \eqn{p > n}). If \eqn{G < p} then RIC may be too restrictive and it would be better to use \code{log(G) < k < 2log(p)}.
#' @param intercept Should an intercept be included in the model? Default is TRUE. Does not count as a step.
#' @param center Should the columns of the design matrix be centered? Default is TRUE.
#' @param normalize Should the design matrix be normalized? Default is TRUE.
#' @param aicstop Early stopping if AIC increases. Default is 0 corresponding to no early stopping. Positive integer values specify the number of times the AIC is allowed to increase in a row, e.g. with \code{aicstop = 2} the algorithm will stop if the AIC criterion increases for 2 steps in a row. The default of \code{\link[stats]{step}} corresponds to \code{aicstop = 1}.
#' @param verbose Print out progress along the way? Default is FALSE.
#' @return An object of class "groupfs" containing information about the sequence of models in the forward stepwise algorithm. Call the function \code{\link{groupfsInf}} on this object to compute selective p-values.
#' @examples
#' x = matrix(rnorm(20*40), nrow=20)
#' index = sort(rep(1:20, 2))
#' y = rnorm(20) + 2 * x[,1] - x[,4]
#' fit = groupfs(x, y, index, maxsteps = 5)
#' out = groupfsInf(fit)
#' out
#' @seealso \code{\link{groupfsInf}}, \code{\link{factorDesign}}.
groupfs <- function(x, y, index, maxsteps, sigma = NULL, k = 2, intercept = TRUE, center = TRUE, normalize = TRUE, aicstop = 0, verbose = FALSE) {

  if (missing(index)) stop("Missing argument: index.")
  p <- ncol(x)
  n <- nrow(x)

  # Group labels
  labels <- unique(index)
  G <- length(labels)
  inactive <- labels
  active <- c()

  if (missing(maxsteps) || maxsteps >= min(n, G)) maxsteps <- min(n-1, G)
  checkargs.xy(x=x, y=y)
  checkargs.groupfs(x, index, maxsteps)
  if (maxsteps > G) stop("maxsteps is larger than number of groups")
  gsizes <- sort(rle(sort(index))$lengths, decreasing = TRUE)
  if (sum(gsizes[1:maxsteps]) >= nrow(x)) {
      maxsteps <- max(which(cumsum(gsizes) < nrow(x)))
      warning(paste("If the largest groups are included the model will be saturated/overdetermined. To prevent this maxsteps has been changed to", maxsteps))
  }

  # Initialize copies of data for loop
  by <- mean(y)
  y.update <- y
  if (intercept) y.update <- y - by
  y.last <- y.update

  # Center and scale design matrix
  xscaled <- scaleGroups(x, index, center, normalize)
  xm <- xscaled$xm
  xs <- xscaled$xs
  x.update <- xscaled$x

  x.begin <- x.update
  y.begin <- y.update
  stopped <- FALSE
  # Store all projections computed along the path
  terms = projections = maxprojs = aicpens = maxpens = cumprojs = vector("list", maxsteps)

  # Store other information from each step
  path.info <- data.frame(imax=integer(maxsteps), df=integer(maxsteps), AIC=numeric(maxsteps), RSS=numeric(maxsteps), RSSdrop=numeric(maxsteps), chisq=numeric(maxsteps))

  modelrank <- as.numeric(intercept)
  if (is.null(sigma)) {
      modelrank <- modelrank + 1
      aic.begin <- aic.last <- n*(log(2*pi) + log(mean(y.update^2))) + k * (n + modelrank)
  } else {
      aic.begin <- aic.last <- sum(y.update^2)/sigma^2 - n + k * modelrank
  }
  if (verbose) print(paste0("Start:  AIC=", round(aic.begin, 3)), quote = FALSE)

  # Begin main loop
  for (step in 1:maxsteps) {

    added <- add1.groupfs(x.update, y.update, index, labels, inactive, k, sigma)

    # Group to be added
    imax <- added$imax
    inactive <- setdiff(inactive, imax)
    active <- union(active, imax)
    inactive.inds <- which(!index %in% active)

    # Rank of group
    modelrank <- modelrank + added$df

    # Stop without adding if model has become saturated
    if (modelrank >= n) {
        stop("Saturated model. Abandon ship!")
    }

    # Regress added group out of y and inactive x
    P.imax <- added$maxproj %*% t(added$maxproj)
    P.imax <- diag(rep(1, n)) - P.imax
    y.update <- P.imax %*% y.update
    x.update[, inactive.inds] <- P.imax %*% x.update[, inactive.inds]

    # Compute AIC
    if (is.null(sigma)) {
        added$AIC <- n * log(added$maxterm/n) - k * added$df + n*log(2*pi) + k * (n + modelrank)
    } else {
        added$AIC <- sum(y.update^2)/sigma^2 - n + k * modelrank
    }

    projections[[step]] <- added$projections
    maxprojs[[step]] <- added$maxproj
    aicpens[[step]] <- added$aicpens
    maxpens[[step]] <- added$maxpen
    if (step == 1) cumprojs[[step]] <- P.imax
    if (step > 1) cumprojs[[step]] <- P.imax %*% cumprojs[[step-1]]
    terms[[step]] <- added$terms

    # Compute RSS for unadjusted chisq p-values
    added$RSS <- sum(y.update^2)
    scale.chisq <- 1

    added$RSSdrop <- sum((y.last - y.update)^2)
    added$chisq <- pchisq(added$RSSdrop/scale.chisq, lower.tail=FALSE, df = added$df)
    y.last <- y.update

    # Projections are stored separately
    step.info <- data.frame(added[-c(3:(length(added)-4))])
    path.info[step, ] <- step.info

    if (verbose) print(round(step.info, 3))

    if (aicstop > 0 && step < maxsteps && step >= aicstop && aic.last < added$AIC) {
        if (all(diff(c(aic.begin, path.info$AIC)[(step+1-aicstop):(step+1)]) > 0)) {

            if (is.null(sigma)) {
                added$AIC <- n * log(added$maxterm/n) - k * added$df + n + n*log(2*pi) + k * modelrank
            } else {
                added$AIC <- sum(y.update^2)/sigma^2 - n + k * modelrank
            }

            path.info <- path.info[1:step, ]
            projections[(step+1):maxsteps] <- NULL
            maxprojs[(step+1):maxsteps] <- NULL
            aicpens[(step+1):maxsteps] <- NULL
            maxpens[(step+1):maxsteps] <- NULL
            cumprojs[(step+1):maxsteps] <- NULL
            terms[(step+1):maxsteps] <- NULL
            maxsteps <- step
            stopped <- TRUE
            break
        }
    }
    aic.last <- added$AIC
  }

  # Is there a better way of doing this?
  # Use some projections already computed?
  beta <- coef(lm(y.begin ~ x.begin[,index %in% path.info$imax]-1))
  names(beta) <- index[index %in% path.info$imax]

  # Create output object
  value <- list(action = path.info$imax, L = path.info$L, AIC = path.info$AIC, projections = projections, maxprojs = maxprojs, aicpens = aicpens, maxpens = maxpens, cumprojs = cumprojs, log = path.info, index = index, y = y.begin, x = x.begin, coefficients = beta, bx = xm, by = by, sx = xs, sigma = sigma, intercept = intercept, call = match.call(), terms = terms)

  class(value) <- "groupfs"
  attr(value, "center") <- center
  attr(value, "normalize") <- normalize
  attr(value, "labels") <- labels
  attr(value, "maxsteps") <- maxsteps
  attr(value, "sigma") <- sigma
  attr(value, "k") <- k
  attr(value, "aicstop") <- aicstop
  attr(value, "stopped") <- stopped
  if (is.null(attr(x, "varnames"))) {
    attr(value, "varnames") <- colnames(x)
  } else {
    attr(value, "varnames") <- attr(x, "varnames")
  }
  return(value)
}

#' Add one group to the model in \code{groupfs}.
#'
#' For internal use by \code{\link{groupfs}}.
#'
#' @param xr Design matrix at current step.
#' @param yr Response vector residual at current step.
#' @param index Group membership indicator of length p.
#' @param labels The unique elements of \code{index}.
#' @param inactive Labels of inactive groups.
#' @param k Multiplier of model size penalty, use \code{k = 2} for AIC, \code{k = log(n)} for BIC, or \code{k = log(p)} for RIC.
#' @param sigma Estimate of error standard deviation for use in AIC criterion. This determines the relative scale between RSS and the degrees of freedom penalty. See \code{\link{extractAIC}} for details.
#' @return Index \code{imax} of added group, value \code{L} of maximized negative AIC, lists of projection matrices defining quadratic model selection event.
add1.groupfs <- function(xr, yr, index, labels, inactive, k, sigma = NULL) {

  # Use characters to avoid issues where
  # list() populates NULL lists in the positions
  # of the active variables
  ### Question for later: does this slow down lapply?
  keys = as.character(inactive)
  n <- nrow(xr)

  # Compute sums of squares to determine which group is added
  # penalized by rank of group if k > 0
  projections = aicpens = terms = vector("list", length(keys))
  names(projections) = names(terms) = names(aicpens) = keys
  for (key in keys) {
      inds <- which(index == key)
      xi <- xr[,inds]
      ui <- svdu_thresh(xi)
      dfi <- ncol(ui)
      projections[[key]] <- ui
      uy <- t(ui) %*% yr
      if (is.null(sigma)) {
          aicpens[[key]] <- exp(k*dfi/n)
          terms[[key]] <- (sum(yr^2) - sum(uy^2)) * aicpens[[key]]
      } else {
          aicpens[[key]] <- sigma^2 * k * dfi
          terms[[key]] <- (sum(yr^2) - sum(uy^2)) + aicpens[[key]]
      }
  }

  # Maximizer = group to be added
  terms.optind <- which.min(terms)
  imax <- inactive[terms.optind]
  optkey <- which(keys == imax)
  maxproj <- projections[[optkey]]
  maxpen <- aicpens[[optkey]]
  maxterm <- terms[[optkey]]
  projections[[optkey]] <- NULL
  aicpens[[optkey]] <- NULL

  return(list(imax=imax, df = ncol(maxproj), projections = projections, maxproj = maxproj, aicpens = aicpens, maxpen = maxpen, maxterm = maxterm, terms = terms))
}

# -----------------------------------------------------------

#' Compute selective p-values for a model fitted by \code{groupfs}.
#'
#' Computes p-values for each group of variables in a model fitted by \code{\link{groupfs}}. These p-values adjust for selection by truncating the usual \eqn{\chi^2} statistics to the regions implied by the model selection event. If the \code{sigma} to \code{\link{groupfs}} was NULL then groupfsInf uses truncated \eqn{F} statistics instead of truncated \eqn{\chi}. The \code{sigma} argument to groupfsInf allows users to override and use \eqn{\chi}, but this is not recommended unless \eqn{\sigma} can be estimated well (i.e. \eqn{n > p}).
#'
#' @param obj Object returned by \code{\link{groupfs}} function
#' @param sigma Estimate of error standard deviation. Default is NULL and in this case groupfsInf uses the value of sigma specified to \code{\link{groupfs}}.
#' @param verbose Print out progress along the way? Default is TRUE.
#' @return An object of class "groupfsInf" containing selective p-values for the fitted model \code{obj}. For comparison with \code{\link{fsInf}}, note that the option \code{type = "active"} is not available.
#'
#' \describe{
#'   \item{vars}{Labels of the active groups in the order they were included.}
#'   \item{pv}{Selective p-values computed from appropriate truncated distributions.}
#'   \item{sigma}{Estimate of error variance used in computing p-values.}
#'   \item{TC or TF}{Observed value of truncated \eqn{\chi} or \eqn{F}.}
#'   \item{df}{Rank of group of variables when it was added to the model.}
#'   \item{support}{List of intervals defining the truncation region of the corresponding statistic.}
#' }
groupfsInf <- function(obj, sigma = NULL, verbose = TRUE) {

  if (!is.null(obj$cvobj) && attr(obj, "stopped")) {
      stop("Cross-validation and early stopping cannot be used simultaneously.")
      # This shouldn't happen in the first place!
      # (it wouldn't anyway unless someone tries to trick it)
  }

  n <- nrow(obj$x)
  p <- ncol(obj$x)
  maxsteps <- attr(obj, "maxsteps")
  k <- attr(obj, "k")
  index <- obj$index
  x <- obj$x
  y <- obj$y
  Ep <- sum(index %in% obj$action)

  pvals = dfs = dfs2 = Tstats = numeric(maxsteps)
  supports <- list()

  if (!is.null(sigma)) {
      type <- "TC"
      if (!is.null(obj$sigma)) {
          cat(paste("Using specified value", sigma, "for sigma in place of the value", obj$sigma, "used by groupfs()\n"))
      }
  } else {
      if (is.null(obj$sigma)) {
          type <- "TF"
          Pf <- svdu_thresh(obj$x[,which(obj$index %in% obj$action), drop = FALSE])
          dffull <- ncol(Pf)
          df2 <- n - dffull - obj$intercept - 1
          Pfull <- Pf %*% t(Pf)
      } else {
          type <- "TC"
          sigma <- obj$sigma
      }
  }

  # Compute p-value for each active group
  for (j in 1:maxsteps) {
    i <- obj$action[j]
    if (verbose) {
        string <- paste0("Step ", j, "/", attr(obj, "maxsteps"), ": computing P-value for group ", i)
        if (!is.null(obj$cvobj)) string <- paste0(string, ", including constraints from cross-validation")
        if (attr(obj, "stopped")) string <- paste0(string, ", including constraints from AICstop")
        cat(paste(string, "\n"))
    }

    if (type == "TC") {
        # Form projection onto active set minus i
        # and project x_i orthogonally
        x_i <- obj$x[,which(obj$index == i), drop = FALSE]
        if (length(obj$action) > 1) {
            minus_i <- setdiff(obj$action, i)
            x_minus_i <- svdu_thresh(obj$x[,which(obj$index %in% minus_i), drop = FALSE])
            x_i <- x_i - x_minus_i %*% t(x_minus_i) %*% x_i
        }

        # Project y onto what remains of x_i
        Ugtilde <- svdu_thresh(x_i)
        R <- t(Ugtilde) %*% obj$y
        TC <- sqrt(sum(R^2))
        eta <- Ugtilde %*% R / TC
        Z <- obj$y - eta * TC
        dfi <- ncol(Ugtilde)
        Tstats[j] <- TC
        dfs[j] <- dfi

        ydecomp <- list(Z=Z, eta=eta)

    } else {

        if (length(obj$action) > 1) {
            minus_i <- setdiff(obj$action, i)
            Psub <- svdu_thresh(obj$x[,which(obj$index %in% minus_i), drop = FALSE])
            Z <- Psub %*% t(Psub) %*% obj$y
            df1 <- dffull - ncol(Psub)
        } else {
            Z <- rep(0, n)
            df1 <- dffull + obj$intercept + 1
        }

        C <- df1/df2
        R1 <- obj$y - Z
        R2 <- obj$y - Pfull %*% obj$y
        R1sq <- sum(R1^2)
        R2sq <- sum(R2^2)
        R <- sqrt(R1sq)
        delta <- R1-R2
        Vdelta <- delta/sqrt(sum(delta^2))
        V2 <- R2/sqrt(R2sq)
        TF <- (R1sq-R2sq)/(C*R2sq)
        Tstats[j] <- TF
        dfs[j] <- df1

        ydecomp <- list(R=R, Z=Z, Vd=Vdelta, V2=V2, C=C)

    }

    intervallist <- truncationRegion(obj, ydecomp, type)

    # Additional constraints from cross-validation?
    if (!is.null(obj$cvobj)) {
        intervallist <- c(intervallist, do.call(c,
                          lapply(obj$cvobj, function(cvf) {
                              if (type == "TC") {
                                  ydecomp <- list(R=R[-cvf$fold], eta=eta[-cvf$fold], Z=Z[-cvf$fold])
                              } else {
                                  ydecomp <- list(R=R, Z=Z[-cvf$fold], Vd=Vdelta[-cvf$fold], V2=V2[-cvf$fold], C=C) # C correct?
                              }
                              truncationRegion(cvf, ydecomp, type)
                          })))
        intervallist <- c(intervallist,
                          lapply(obj$cvquad, function(cvquad) {
                              if (type == "TC") {
                                  etacvquad <- t(eta) %*% cvquad
                                  A <- etacvquad %*% eta
                                  B <- 2 * etacvquad %*% Z
                                  C <- t(Z) %*% cvquad %*% Z
                                  quadratic_roots(A, B, C, tol = 1e-15)
                              } else {

                                  zcvquad <- t(Z) %*% cvquad
                                  vdcvquad <- t(Vdelta) %*% cvquad
                                  v2cvquad <- t(V2) %*% cvquad
                                  x0 <- zcvquad %*% Z
                                  x1 <- 2*R*zcvquad %*% Vdelta
                                  x2 <- 2*R*zcvquad %*% V2
                                  x12 <- 2*R^2*vdcvquad %*% V2
                                  x11 <- R^2*vdcvquad %*% Vdelta
                                  x22 <- R^2*v2cvquad %*% V2
                                  TF_roots(R, C, coeffs = list(x0=x0, x1=x1, x2=x2, x12=x12, x11=x11, x22=x22))
                              }
                          }))
    }

    # Additional constraints from AIC stopping
    if (attr(obj, "stopped")) {
        aicintervals <- vector("list", maxsteps)
        aicstop <- attr(obj, "aicstop")
        if (type == "TC") {
            pen0 <- k * obj$intercept
            aic.begin <- aic.last <- sum(obj$y^2)/sigma^2 - n + k * obj$intercept
        } else {
            pen0 <- exp(k * (1+obj$intercept)/n)
            aic.begin <- n*(log(2*pi) + log(mean(obj$y^2))) + k * (1 + n + obj$intercept)
        }
        AICs <- c(aic.begin, obj$AIC)

        ulist <- c(list(matrix(0, n, 1)), obj$maxprojs)
        penlist <- c(pen0, obj$maxpens)
        zlist <- vector("list", maxsteps+1)
        zlist[[1]] <- zlist[[2]] <- Z
        if (type == "TC") {
            etalist <- vector("list", maxsteps+1)
            etalist[[1]] <- etalist[[2]] <- eta
        } else {
            vdlist <- v2list <- vector("list", maxsteps+1)
            vdlist[[1]] <- vdlist[[2]] <- Vdelta
            v2list[[1]] <- v2list[[2]] <- V2
        }
        if (maxsteps > 1) {
            for (step in 1:(maxsteps-1)) {
                cproj <- obj$cumprojs[[step]]
                zlist[[step+2]] <- cproj %*% Z
                if (type == "TC") {
                    etalist[[step+2]] <- cproj %*% eta
                } else {
                    vdlist[[step+2]] <- cproj %*% Vdelta
                    v2list[[step+2]] <- cproj %*% V2
                }
            }
        }

        for (step in 1:maxsteps) {
            # Compare AIC at s+1 to AIC at s
            # roots() functions assume g indexes smaller AIC
            # this is step+1 until the last step
            peng <- penlist[[step+1]]
            Ug <- ulist[[step+1]]
            Uh <- ulist[[step]]
            Zg <- zlist[[step+1]]
            Zh <- zlist[[step]]

            if (type == "TC") {
                penh <- 0
                etag <- etalist[[step+1]]
                etah <- etalist[[step]]
                coeffs <- quadratic_coefficients(obj$sigma, Ug, Uh, peng, penh, etag, etah, Zg, Zh)

                if (AICs[step] < AICs[step+1]) {
                    coeffs <- lapply(coeffs, function(coeff) -coeff)
                }

                intstep <- quadratic_roots(coeffs$A, coeffs$B, coeffs$C, tol = 1e-15)

            } else {
                penh <- 1
                Vdg <- vdlist[[step+1]]
                Vdh <- vdlist[[step]]
                V2g <- v2list[[step+1]]
                V2h <- v2list[[step]]
                coeffs <- TF_coefficients(R, Ug, Uh, peng, penh, Zg, Zh, Vdg, Vdh, V2g, V2h)

                if (AICs[step] < AICs[step+1]) {
                    coeffs <- lapply(coeffs, function(coeff) -coeff)
                }

                intstep <- TF_roots(R, C, coeffs)
            }

            aicintervals[[step]] <- intstep
        }
        intervallist <- c(intervallist, aicintervals)
    }

    # Compute intersection:
    region <- do.call(interval_union, intervallist)
    region <- interval_union(region, Intervals(c(-Inf,0)))
    E <- interval_complement(region, check_valid = FALSE)

    if (length(E) == 0) {
        stop(paste("Empty support at step", j))
    }
    supports[[j]] <- E

    # E is now potentially a union of intervals
    if (type == "TC") {
        pvals[j] <- TC_surv(TC, sigma, dfi, E)
    } else {
        # write TF_surv function first
        pvals[j] <- TF_surv(TF, df1, df2, E)
    }

  }

  if (any(is.nan(pvals))) {
      nanp <- which(is.nan(pvals))
      pvals[nanp] <- 0
      warning(paste0("P-value NaNs of the form 0/0 converted to 0 for group(s) ", paste(obj$action[nanp], collapse=","), ". This typically occurs for numerical reasons in the presence of a large signal-to-noise ratio."))
  }

  names(pvals) <- obj$action
  out <- list(vars = obj$action, pv=pvals)
  if (type == "TC") {
      out$TC <- Tstats
      out$sigma <- sigma
  } else {
      out$TF <- Tstats
      out$df2 <- df2
  }
  out$df <- dfs
  out$support <- supports
  class(out) <- "groupfsInf"
  if (!is.null(attr(obj, "varnames"))) {
      attr(out, "varnames") <- attr(obj, "varnames")
  }
  return(out)
}

# -----------------------------------------------------------

TC_surv <- function(TC, sigma, df, E) {
    if (length(E) == 0) {
        stop("Empty TC support")
    }

    # Sum truncated cdf over each part of E
    denom <- do.call(sum, lapply(1:nrow(E), function(v) {
      tchi_interval(E[v,1], E[v,2], sigma, df)
    }))

    # Sum truncated cdf from observed value to max of
    # truncation region
    numer <- do.call(sum, lapply(1:nrow(E), function(v) {
      lower <- E[v,1]
      upper <- E[v,2]
      if (upper > TC) {
        # Observed value is left of this interval's right endpoint
        if (lower < TC) {
          # Observed value is in this interval
          return(tchi_interval(TC, upper, sigma, df))
        } else {
          # Observed value is not in this interval
          return(tchi_interval(lower, upper, sigma, df))
        }
      } else {
        # Observed value is right of this entire interval
        return(0)
      }
    }))

    # Survival function
    value <- numer/denom
    # Force p-value to lie in the [0,1] interval
    # in case of numerical issues
    value <- max(0, min(1, value))
    value
}

tchi_interval <- function(lower, upper, sigma, df) {
  a <- (lower/sigma)^2
  b <- (upper/sigma)^2
  if (b == Inf) {
      integral <- pchisq(a, df, lower.tail = FALSE)
  } else {
      integral <- pchisq(b, df) - pchisq(a, df)
  }
  if ((integral < .Machine$double.eps) && (b < Inf)) {
      integral <- num_int_chi(a, b, df)
  }
  return(integral)
}

num_int_chi <- function(a, b, df, nsamp = 10000) {
  grid <- seq(from=a, to=b, length.out=nsamp)
  integrand <- dchisq(grid, df)
  return((b-a)*mean(integrand))
}

TF_surv <- function(TF, df1, df2, E) {
    if (length(E) == 0) {
        stop("Empty TF support")
    }

    # Sum truncated cdf over each part of E
    denom <- do.call(sum, lapply(1:nrow(E), function(v) {
      TF_interval(E[v,1], E[v,2], df1, df2)
    }))

    # Sum truncated cdf from observed value to max of
    # truncation region
    numer <- do.call(sum, lapply(1:nrow(E), function(v) {
      lower <- E[v,1]
      upper <- E[v,2]
      if (upper > TF) {
        # Observed value is left of this interval's right endpoint
        if (lower < TF) {
          # Observed value is in this interval
          return(TF_interval(TF, upper, df1, df2))
        } else {
          # Observed value is not in this interval
          return(TF_interval(lower, upper, df1, df2))
        }
      } else {
        # Observed value is right of this entire interval
        return(0)
      }
    }))

    # Survival function
    value <- numer/denom
    # Force p-value to lie in the [0,1] interval
    # in case of numerical issues
    #value <- max(0, min(1, value))
    value
}

TF_interval <- function(lower, upper, df1, df2) {
  a <- lower
  b <- upper
  if (b == Inf) {
      integral <- pf(a, df1, df2, lower.tail = FALSE)
  } else {
      integral <- pf(b, df1, df2) - pf(a, df1, df2)
  }
  if ((integral < .Machine$double.eps) && (b < Inf)) {
      integral <- num_int_F(a, b, df1, df2)
  }
  return(integral)
}

num_int_F <- function(a, b, df1, df2, nsamp = 10000) {
  grid <- seq(from=a, to=b, length.out=nsamp)
  integrand <- df(grid, df1, df2)
  return((b-a)*mean(integrand))
}

#' Center and scale design matrix by groups
#'
#' For internal use by \code{\link{groupfs}}.
#'
#' @param x Design matrix.
#' @param index Group membership indicator of length p.
#' @param center Center groups, default is TRUE.
#' @param normalize Scale groups by Frobenius norm, default is TRUE.
#' @return
#' \describe{
#'   \item{x}{Optionally centered/scaled design matrix.}
#'   \item{xm}{Means of groups in original design matrix.}
#'   \item{xs}{Frobenius norms of groups in original design matrix.}
#' }
scaleGroups <- function(x, index, center = TRUE, normalize = TRUE) {
  keys <- unique(index)
  xm <- rep(0, ncol(x))
  xs <- rep(1, ncol(x))

  for (j in keys) {
    inds <- which(index == j)
    if (center) {
        xmj <- mean(x[, inds])
        xm[inds] <- xmj
        x[, inds] <- x[, inds] - xmj
    }
    normsq <- sum(x[, inds]^2)
    xsj <- sqrt(normsq)
    xs[inds] <- xsj
    if (xsj > 0) {
        if (normalize) x[, inds] <- x[, inds] / xsj
    } else {
        stop(paste("Design matrix contains identically zero group of variables:", j))
    }
  }
  return(list(x=x, xm=xm, xs=xs))
}

#' Expand a data frame with factors to form a design matrix with the full binary encoding of each factor.
#'
#' When using \code{\link{groupfs}} with factor variables call this function first to create a design matrix.
#'
#' @param df Data frame containing some columns which are \code{factors}.
#' @return List containing
#' \describe{
#'   \item{x}{Design matrix, the first columns contain any numeric variables from the original date frame.}
#'   \item{index}{Group membership indicator for expanded matrix.}
#' }
#' @examples
#' \dontrun{
#' fd = factorDesign(warpbreaks)
#' y = rnorm(nrow(fd$x))
#' fit = groupfs(fd$x, y, fd$index, maxsteps=2, intercept=F)
#' pvals = groupfsInf(fit)
#' }
factorDesign <- function(df) {
    factor.inds <- sapply(df[1,], is.factor)
    factor.labels <- which(factor.inds)
    nfacs <- sum(factor.inds)
    nlevs <- sapply(df[1,factor.inds], function(fac) nlevels(fac))
    totnlevs <- sum(nlevs)
    num.num = indcounter = ncol(df) - nfacs
    x <- matrix(nrow=nrow(df), ncol = totnlevs + num.num)
    colnames(x) <- 1:ncol(x)
    index <- integer(ncol(x))
    varnames <- character(ncol(df))
    if (num.num > 0) {
        x[,1:num.num] <- df[, !factor.inds]
        varnames[1:num.num] = colnames(x)[1:num.num] <- colnames(df)[1:num.num]
        index[1:num.num] <- 1:num.num
        indcounter <- indcounter + num.num - 1
    }
    for (j in 1:nfacs) {
        submat <- model.matrix(~ df[, factor.labels[j]] - 1)
        indcounter <- indcounter+1
        submatinds <- indcounter:(indcounter+nlevs[j]-1)
        indcounter <- indcounter + nlevs[j] - 1
        colnames(x)[submatinds] <- paste0(colnames(df)[num.num + j], ":", 1:nlevs[j])
        varnames[num.num + j] <- colnames(df)[num.num + j]
        x[,submatinds] <- submat
        index[submatinds] <- num.num + j
    }
    attr(x, "varnames") <- varnames
    return(list(x = x, index = index))
}

svdu_thresh <- function(x) {
    svdx <- svd(x)
    inds <- svdx$d > svdx$d[1] * sqrt(.Machine$double.eps)
    return(svdx$u[, inds, drop = FALSE])
}

flatten <- function(L) {
    if (is.list(L[[1]])) return(unlist(L, recursive=FALSE))
    return(L)
}

print.groupfs <- function(x, ...) {
    cat("\nSequence of added groups:\n")
    nsteps = length(x$action)
    action <- x$action
    vnames <- attr(x, "varnames")
    if (length(vnames) > 0) action <- vnames[action]
    tab = data.frame(Group = action, Rank = x$log$df, RSS = round(x$log$RSS, 3), AIC = round(x$log$AIC, 3))
    rownames(tab) = 1:nsteps
    print(tab)
    cat("\nUse groupfsInf() to compute P-values\n")
    invisible()
}


coef.groupfs <- function(object, ...) {
    return(object$coefficients)
}

#' @name predict.groupfs
#' @aliases predict.groupfs
#' @aliases coef.groupfs
#'
#' @title Prediction and coefficient functions for \code{\link{groupfs}}.
#'
#' Make predictions or extract coefficients from a groupfs forward stepwise object.
#'
#' @param object Object returned by a call to \code{\link{groupfs}}.
#' @param newx Matrix of x values at which the predictions are desired. If NULL, the x values from groupfs fitting are used.
#' @return A vector of predictions or a vector of coefficients.
predict.groupfs <- function(object, newx) {
    beta <- coef.groupfs(object)
    if (missing(newx)) {
        newx = object$x
    } else {
        newx <- scaleGroups(newx, object$index, attr(object, "center"), attr(object, "normalize"))$x
    }
    return(newx[, object$index %in% object$action] %*% beta + ifelse(object$intercept, object$by, 0))
}

print.groupfsInf <- function(x, ...) {
    if (!is.null(x$sigma)) {
        isTF <- FALSE
        Tstat <- x$TC
        cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n", x$sigma))
    } else {
        isTF <- TRUE
        Tstat <- x$TF
    }
    action <- x$vars
    vnames <- attr(x, "varnames")
    if (length(vnames) > 0) action <- vnames[action]
    tab = data.frame(Group = action, Pvalue = round(x$pv, 3),
        TC = round(Tstat, 3),
        df = x$df, Size = round(unlist(lapply(lapply(x$support, size), sum)), 3),
        Ints = unlist(lapply(x$support, nrow)), Min =round(unlist(lapply(x$support, min)), 3),
        Max = round(unlist(lapply(x$support, max)), 3))
    rownames(tab) = 1:length(x$vars)
    if (isTF) names(tab)[3] <- "TF"
    print(tab)
    cat("\nInts is the number of intervals in the truncated chi selection region and Size is the sum of their lengths. Min and Max are the lowest and highest endpoints of the truncation region. No confidence intervals are reported by groupfsInf.\n")
    invisible()
}

checkargs.groupfs <- function(x, index, maxsteps) {
    if (length(index) != ncol(x)) stop("Length of index does not match number of columns of x")
    if ((round(maxsteps) != maxsteps) || (maxsteps <= 0)) stop("maxsteps must be an integer > 0")
}
