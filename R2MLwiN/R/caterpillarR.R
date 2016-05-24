#' Draws caterpillar plots of the residuals from a mlwinfitIGLS-class or 
#' mlwinfitMCMC-class object, at a chosen level of the multilevel model.
#' 
#' Uses \code{\link[lattice]{qqmath}} in the \pkg{lattice} package to draw Quantile-Quantile plots
#' of the residuals at a chosen level of a multilevel model against a
#' theoretical distribution.
#'  
#' @param resi A \code{\link{mlwinfitIGLS-class}} or 
#' \code{\link{mlwinfitMCMC-class}} object (model must be fitted specifying
#' \code{resi.store = TRUE} and including \code{'variances'} or \code{'sampling'} in the list of
#' \code{resioptions} (both in \code{estoptions}) to not return an error).
#' @param lev An integer scalar specifying the level of a multilevel model
#' for which to produce a plot for.
#' 
#' @return See \code{\link[lattice]{qqmath}}.
#' 
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#' 
#' @seealso \code{\link{caterpillar}}, \code{\link[lattice]{qqmath}}
#' 
#' @examples
#' 
#' \dontrun{
#' library(R2MLwiN)
#' # NOTE: if MLwiN not saved where R2MLwiN defaults to:
#' # options(MLwiN_path = 'path/to/MLwiN vX.XX/')
#' # If using R2MLwiN via WINE, the path may look like:
#' # options(MLwiN_path = '/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN vX.XX/') 
#' 
#' # Example using tutorial dataset
#' data(tutorial, package = 'R2MLwiN')
#' mymodel <- runMLwiN(normexam ~ 1 + (1 | school) + (1 | student),
#'                     estoptions = list(resi.store = TRUE),
#'                     data = tutorial)
#' # Caterpillar plot
#' caterpillarR(mymodel['residual'], lev = 2)
#' }
#' 
#' @export
caterpillarR <- function(resi, lev = 2) {
  ## produce caterpillar plots for the random effects(level>=2) using qqmath() method in lattice package Only work
  ## with full covariance specified
  qqmath.ranef.mer <- function(x, data, main = TRUE, ...) {
    prepanel.ci <- function(x, y, se, subscripts, ...) {
      x <- as.numeric(x)
      se <- as.numeric(se[subscripts])
      hw <- 1.96 * se
      list(xlim = range(x - hw, x + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16, ...) {
      panel.grid(h = -1, v = -1)
      panel.abline(v = 0)
      x <- as.numeric(x)
      y <- as.numeric(y)
      se <- as.numeric(se[subscripts])
      panel.segments(x - 1.96 * se, y, x + 1.96 * se, y, col = "black")
      panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(nx) {
      xt <- x[[nx]]
      mtit <- if (main) 
        nx
      if (!is.null(pv <- attr(xt, "postVar"))) {
        d <- dim(pv)
        se <- vapply(seq_len(d[1]), function(i) sqrt(pv[i, i, ]), numeric(d[3]))
        nr <- nrow(xt)
        nc <- ncol(xt)
        ord <- unlist(lapply(xt, order)) + rep((0:(nc - 1)) * nr, each = nr)
        rr <- 1:nr
        ind <- gl(nc, nr, labels = names(xt))
        xyplot(rep(qnorm((rr - 0.5)/nr), nc) ~ unlist(xt)[ord] | ind[ord], se = se[ord], prepanel = prepanel.ci, 
               panel = panel.ci, scales = list(x = list(relation = "free")), ylab = "Standard normal quantiles", 
               xlab = NULL, main = mtit, ...)
      } else {
        qqmath(~values | ind, stack(xt), scales = list(y = list(relation = "free")), xlab = "Standard normal quantiles", 
               ylab = NULL, main = mtit, ...)
      }
    }
    sapply(names(x), f, simplify = FALSE)
  }
  
  if (class(resi) == "mlwinfitIGLS" || class(resi) == "mlwinfitMCMC") {
    myresi <- resi@residual
    if (is.null(myresi)) {
      stop("To generate a caterpillar plot the model must be run with the resi.store option set to TRUE")
    }
  } else {
    if (is.list(resi)) {
      myresi <- resi
    } else {
      if (is.data.frame(resi)) {
        myresi <- as.list(resi)
      } else {
        stop("Invalid resi option specified")
      }
    }
  }
  
  est.names <- names(myresi)[grep(paste("lev_", lev, "_resi_est", sep = ""), names(myresi))]
  if (length(est.names) == 1) {
    est <- as.matrix(na.omit(myresi[[est.names]]))
    colnames(est) <- sub("_resi_est", "", est.names)
    var <- na.omit(myresi[[grep(paste("lev_", lev, "_resi_(var|variance)_", sep = ""), names(myresi))[1]]])
    d1 <- length(est)
    tt <- array(, c(1, 1, d1))
    tt[1, 1, ] <- var
  } else {
    est <- NULL
    for (i in 1:length(est.names)) {
      est <- cbind(est, myresi[[est.names[i]]])
    }
    est <- na.omit(est)
    colnames(est) <- sub("_resi_est", "", est.names)
    tempnames <- sub(paste("lev_", lev, "_resi_est_", sep = ""), "", est.names)
    d1 <- dim(est)[1]
    d2 <- dim(est)[2]
    m <- (d2 * (d2 + 1))/2
    cov.lower <- matrix(NA, m, d1)
    ccount <- 1
    for (i in 1:length(est.names)) {
      for (j in 1:i) {
        if (i == j) {
          tmatch <- grep(paste("lev_", lev, "_resi_(var|variance)_", tempnames[i], sep = ""), names(myresi))[1]
          if (length(tmatch) != 0) {
            cov.lower[ccount, ] <- na.omit(myresi[[tmatch]])
          } else {
            cov.lower[ccount, ] <- rep(0, d1)
          }
        } else {
          tmatch <- grep(paste("lev_", lev, "_resi_cov_", tempnames[i], "_", tempnames[j], sep = ""), names(myresi))
          if (length(tmatch) != 0) {
            cov.lower[ccount, ] <- na.omit(myresi[[tmatch]])
          } else {
            cov.lower[ccount, ] <- rep(0, d1)
          }
        }
        ccount <- ccount + 1
      }
    }
    tt <- array(, c(d2, d2, d1))
    for (x in 1:d1) {
      tt[, , x][upper.tri(tt[, , x], diag = TRUE)] <- cov.lower[, x]
      tt[, , x][lower.tri(tt[, , x], diag = FALSE)] <- t(tt[, , x])[lower.tri(tt[, , x], diag = FALSE)]
    }
  }
  
  rr <- NULL
  rr$Subject <- data.frame(est)
  attr(rr$Subject, "postVar") <- tt
  class(rr) <- "ranef.mer"
  qqmath.ranef.mer(rr)
} 
