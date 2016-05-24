#' Draws predicted curves (lines) using estimates from the fixed part of a
#' fitted model.
#'
#' This function draws predicted curves (lines) against an explanatory variable
#' for each category of a categorical variable.
#'
#' @param object Either an \code{\link{mlwinfitIGLS-class}} or \code{\link{mlwinfitMCMC-class}} object.
#' @param indata A data.frame object containing the data. If not specified, data is extracted from
#' the \code{object}.
#' @param xname The name of variable to be plotted.
#' @param group A character string or a sequence of length equivalent to rows of data to plot.
#' \code{group = NULL} by default.
#' @param legend A logical value indicating whether a legend for \code{group}
#' is to be added.
#' @param legend.space A character string specifies one of the four sides,
#' which can be one of \code{'top'}, \code{'bottom'}, \code{'left'} and \code{'right'}. Default,
#' \code{legend.space = 'top'}.
#' @param legend.ncol An integer specifies a number of columns, possibly
#' divided into blocks, each containing some rows. Default,
#' \code{legend.ncol = 2}.
#' @param ...  Other arguments to be passed to \code{\link[lattice]{xyplot}}.
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso \code{\link{predLines}}
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
#' ## Read alevchem data
#' data(alevchem, package = "R2MLwiN")
#'
#' alevchem$gcseav <- double2singlePrecision(alevchem$gcse_tot/alevchem$gcse_no - 6)
#' # Avoids warning when fitting factor as continuous response:
#' alevchem$a_point_num <- as.numeric(alevchem$a_point)
#'
#' ## Example: A-level Chemistry
#' (mymodel <- runMLwiN(a_point_num ~ 1 + gcseav + I(gcseav^2) + I(gcseav^3)
#'                      + gender + (1 | pupil), estoptions = list(EstM = 1,  resi.store = TRUE),
#'                      data = alevchem))
#'                      
#' predCurves(mymodel, xname = "gcseav", group = "genderfemale")
#' }
#'
#' @export
predCurves <- function(object, indata = NULL, xname, group = NULL, legend = TRUE, legend.space = "top", legend.ncol = 2, 
                       ...) {
  ## This function is to draw predicted lines using fixed part estimates
  
  cls <- class(object)
  if (!cls %in% c("mlwinfitIGLS", "mlwinfitMCMC")) 
    stop("need a \"mlwinfitIGLS\" or \"mlwinfitMCMC\" class object")
  
  FP <- object@FP
  if (is.null(indata)) {
    indata <- object@data
  }
  if (!xname %in% colnames(indata)) {
    stop(paste(xname, " does not exist in the data"))
  }
  if (!is.null(group)) {
    if (is.character(group)) 
      group <- indata[[group]]
    if (!is.factor(group)) 
      group <- as.factor(group)
  }
  
  fp.names <- sub("FP_", "", names(FP))
  tval <- 0
  for (i in 1:length(fp.names)) {
    if (is.factor(indata[[fp.names[i]]])) {
      indata[[fp.names[i]]] <- as.integer(indata[[fp.names[i]]]) - 1
    }
    tval <- tval + as.numeric(indata[[fp.names[i]]]) * FP[i]
  }
  
  pred.min <- min(tval)
  pred.max <- max(tval)
  x <- indata[[xname]]
  x.min <- min(x)
  x.max <- max(x)
  
  if (legend && length(group)) {
    key <- list(lines = Rows(trellis.par.get("superpose.line"), 1:nlevels(group)), text = list(lab = levels(group)), 
                space = legend.space, columns = legend.ncol)
  } else {
    key <- NULL
  }
  
  if (!is.null(group)) {
    levs <- levels(group)
    nlev <- length(levs)
    trellis.obj <- xyplot(tval ~ x, prepanel = function(x, y, ...) {
      list(xlim = c(x.min, x.max), ylim = c(pred.min, pred.max))
    }, groups = group, panel = function(x, y, groups, ...) {
      col <- Rows(trellis.par.get("superpose.line"), 1:nlev)$col
      for (i in 1:nlev) {
        ypred <- y[groups == levs[i]]
        panel.xyplot(x = sort(x[groups == levs[i]]), y = ypred[order(x[groups == levs[i]])], col = col[i], type = "l", ...)
      }
    }, key = key, ylab = "ypred", xlab = xname, ...)
  } else {
    trellis.obj <- xyplot(tval ~ x, prepanel = function(x, y, ...) {
      list(xlim = c(x.min, x.max), ylim = c(pred.min, pred.max))
    }, panel = function(x, y, ...) {
      panel.xyplot(x = sort(x), y = y[order(x)], type = "l", ...)
    }, ylab = "ypred", xlab = xname, ...)
  }
  print(trellis.obj)
  invisible(trellis.obj)
} 
