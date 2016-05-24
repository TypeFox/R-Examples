#' Draws predicted lines using a fitted model object
#'
#' This function draws predicted lines against an explanatory variable for
#' selected groups at a higher (>=2) level.
#'
#' @param object Either an \code{\link{mlwinfitIGLS-class}} or \code{\link{mlwinfitMCMC-class}}
#' object.
#' @param indata A data.frame object containing the data. If not specified, data is extracted from
#' the \code{object}.
#' @param xname The name of the variable to be plotted.
#' @param lev A digit indicating the level (of the multilevel model) at which
#' to plot.
#' @param selected A vector specifying groups to selectively plot at the level
#' specified in \code{lev}. If \code{selected = NULL}, then all groups at that
#' level are included.
#' @param probs A numeric vector of probabilities with values in \code{[0, 1]}
#' used to calculate the lower and upper quantiles from which the error bars
#' are plotted. Currently, this is only available for an \code{\link{mlwinfitMCMC-class}} object.
#' @param legend A logical value indicating whether a legend is to be added.
#' @param legend.space A character string specifies one of the four sides,
#' which can be one of \code{'top'}, \code{'bottom'}, \code{'left'} and \code{'right'}.  Default,
#' \code{legend.space = 'top'}.
#' @param legend.ncol An integer specifies a number of columns, possibly
#' divided into blocks, each containing some rows. Default,
#' \code{legend.ncol = 2}.
#' @param ...  Other arguments to be pased to \code{\link[lattice]{xyplot}}.
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso \code{\link{predCurves}}
#'
#' @examples
#'
#' \dontrun{
#' library(R2MLwiN)
#' # NOTE: if MLwiN not saved in location R2MLwiN defaults to, specify path via:
#' # options(MLwiN_path = 'path/to/MLwiN vX.XX/')
#' # If using R2MLwiN via WINE, the path may look like this:
#' # options(MLwiN_path = '/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN vX.XX/')
#'
#' ## Example: tutorial
#' data(tutorial, package = "R2MLwiN")
#' (mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student),
#'                      estoptions = list(EstM = 1, resi.store.levs = 2), data = tutorial))
#'
#' predLines(mymodel, xname = "standlrt", lev = 2, selected = c(30, 44, 53, 59),
#'           probs = c(0.025, 0.975))
#' }
#'
#' @export
predLines <- function(object, indata = NULL, xname, lev = 2, selected = NULL, probs = c(0.025, 0.975), legend = TRUE, 
                      legend.space = "top", legend.ncol = 4, ...) {
  
  ## This function is to draw predicted lines at higher levels (level>=2)
  if (lev < 2) 
    stop("lev has to be greater than 1")
  cls <- class(object)
  if (!cls %in% c("mlwinfitIGLS", "mlwinfitMCMC")) 
    stop("need a \"mlwinfitIGLS\" or \"mlwinfitMCMC\" class object")
  if (is.null(indata)) {
    indata <- object[["data"]]
  }
  if (cls == "mlwinfitIGLS") {
    FP <- object@FP
    myresi <- object@residual
    levID <- object@levID
    if (length(myresi) == 0) {
      stop("Residuals were not stored")
    }
    
    categrv <- as.factor(indata[[rev(levID)[lev]]])
    levels(categrv) <- 1:length(levels(categrv))
    categrv <- as.integer(categrv)
    if (is.null(selected)) {
      selected <- unique(categrv)
    }
    
    est.names <- names(myresi)[grep(paste("lev_", lev, "_resi_est", sep = ""), names(myresi))]
    if (length(est.names) == 0) {
      stop("Residuals were not stored at the requested level")
    }
    if (length(est.names) == 1) {
      est0 <- na.omit(myresi[[est.names]])
      if (length(est0) == length(unique(categrv))) {
        est <- as.matrix(est0[categrv])
        colnames(est) <- sub("_resi_est", "", est.names)
      } else {
        stop("The number of groups do not match the number of residual estimates.")
      }
    } else {
      est0 <- NULL
      for (i in 1:length(est.names)) {
        est0 <- cbind(est0, myresi[[est.names[i]]])
      }
      if (nrow(est0) == length(unique(categrv))) {
        est <- as.matrix(est0[categrv, ])
        colnames(est) <- sub("_resi_est", "", est.names)
      } else {
        stop("The number of groups do not match the number of residual estimates.")
      }
    }
    
    rpx.names <- sub(paste("lev_", lev, "_", sep = ""), "", colnames(est))
    fp.names <- sub("FP_", "", names(FP))
    tval <- 0
    for (i in 1:length(fp.names)) {
      if (is.factor(indata[[fp.names[i]]])) {
        indata[[fp.names[i]]] <- as.integer(indata[[fp.names[i]]]) - 1
      }
      tval <- tval + as.numeric(indata[[fp.names[i]]]) * FP[i]
    }
    for (i in 1:length(rpx.names)) {
      if (is.factor(indata[[rpx.names[i]]])) {
        indata[[rpx.names[i]]] <- as.integer(indata[[rpx.names[i]]]) - 1
      }
      tval <- tval + indata[[rpx.names[i]]] * est[, i]
    }
    
    pred.min <- min(tval)
    pred.max <- max(tval)
    pred.diff <- pred.max - pred.min
    x <- indata[[xname]]
    x.min <- min(x)
    x.max <- max(x)
    
    if (legend) {
      key <- list(lines = Rows(trellis.par.get("superpose.line"), 1:length(selected)), text = list(lab = as.character(selected)), 
                  space = legend.space, columns = legend.ncol)
    } else {
      key <- NULL
    }
    
    trellis.obj <- xyplot(tval ~ x, prepanel = function(x, y, ...) {
      list(xlim = c(x.min, x.max), ylim = c(pred.min, pred.max))
    }, groups = categrv, panel = function(x, y, groups, ...) {
      col <- Rows(trellis.par.get("superpose.line"), 1:length(selected))$col
      j <- 1
      for (i in selected) {
        ypred <- y[which(groups == i)]
        panel.xyplot(x = sort(x[which(groups == i)]), y = ypred[order(x[which(groups == i)])], col = col[j], 
                     type = "l", ...)
        j <- j + 1
      }
    }, key = key, ylab = "ypred", xlab = xname, ...)
    print(trellis.obj)
  }
  if (cls == "mlwinfitMCMC") {
    
    ## This function is to draw predicted lines (medians, lower quantiles and upper quantiles) at higher levels
    ## (level>=2)
    resi.chains <- object@resi.chains
    chains <- object@chains
    levID <- object@levID
    if (is.null(resi.chains)) {
      stop("Residual chains were not stored")
    }
    categrv <- indata[[rev(levID)[lev]]]
    if (is.null(selected)) {
      selected <- unique(categrv)
    }
    
    rpx.names <- sub(paste("RP", lev, "_var_", sep = ""), "", colnames(chains)[grep(paste("RP", lev, "_var_", 
                                                                                          sep = ""), colnames(chains))])
    lenrpx <- length(rpx.names)
    if (length(rpx.names) == 0) {
      stop("Residual chains were not stored at the requested level")
    }
    resi.chains <- resi.chains[[paste0("resi_lev", lev)]]
    
    FP.pos <- grep("FP_", colnames(chains))
    fp.names <- sub("FP_", "", colnames(chains)[FP.pos])
    tval <- matrix(0, nrow(indata), nrow(chains))
    
    for (j in selected) {
      for (i in 1:length(fp.names)) {
        fpxdata <- indata[categrv == j, fp.names[i]]
        if (is.factor(fpxdata)) {
          fpxdata <- as.integer(fpxdata) - 1
        }
        tval[categrv == j, ] <- tval[categrv == j, ] + fpxdata %o% chains[, FP.pos[i]]
      }
      for (i in 1:length(rpx.names)) {
        rpxdata <- indata[categrv == j, rpx.names[i]]
        if (is.factor(rpxdata)) {
          rpxdata <- as.integer(rpxdata) - 1
        }
        tval[categrv == j, ] <- tval[categrv == j, ] + rpxdata %o% resi.chains[, paste0("u_", i-1, "_", j)]
      }
    }
    
    quants <- apply(tval, 1, function(x) quantile(x, c(probs[1], 0.5, probs[2])))
    tval.med <- quants[2,]
    tval.low <- quants[1,]
    tval.up <- quants[3,]

    pred.min <- min(tval.low)
    pred.max <- max(tval.up)
    pred.diff <- pred.max - pred.min
    x <- indata[[xname]]
    x.min <- min(x)
    x.max <- max(x)
    
    if (legend) {
      key <- list(lines = Rows(trellis.par.get("superpose.line"), 1:length(selected)), text = list(lab = as.character(selected)), 
                  space = legend.space, columns = legend.ncol)
    } else {
      key <- NULL
    }
    
    trellis.obj <- xyplot(tval.med ~ x, prepanel = function(x, y, ...) {
      list(xlim = c(x.min, x.max), ylim = c(pred.min, pred.max))
    }, groups = categrv, panel = function(x, y, groups, ...) {
      col <- Rows(trellis.par.get("superpose.line"), 1:length(selected))$col
      j <- 1
      for (i in selected) {
        xg <- x[which(groups == i)]
        ypred <- y[which(groups == i)]
        ypred.low <- tval.low[which(groups == i)]
        ypred.up <- tval.up[which(groups == i)]
        panel.xyplot(x = sort(xg), y = ypred[order(xg)], col = col[j], type = "l", ...)
        panel.xyplot(x = sort(xg), y = ypred.low[order(xg)], col = col[j], type = "l", lty = 3, ...)
        panel.xyplot(x = sort(xg), y = ypred.up[order(xg)], col = col[j], type = "l", lty = 3, ...)
        j <- j + 1
      }
    }, key = key, ylab = "ypred", xlab = xname, ...)
    print(trellis.obj)
  }
  invisible(trellis.obj)
} 
