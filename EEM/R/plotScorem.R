#' Plot score matrix for prcomp result based on group
#' 
#' Plot score matrix for \code{\link[stats]{prcomp}} (PCA) result based on group
#' 
#' @inheritParams plotScore
#' @param ncomp maximum number of PC score to plot
#' @param legendtitle legend title
#' @param group a vector of numeric, character or factor class separating 
#' the samples into groups.
#' @param ... additional arguments to be passed on to \code{\link[graphics]{pairs}} 
#' 
#' @return A figure is returned on the graphic device
#' 
#' @seealso
#' \code{\link{pairs}}, \code{\link{plotScore}}
#' 
#' @examples
#' data(applejuice)
#' # country of apple production
#' country <- sapply(strsplit(names(applejuice), split = "-"), "[", 1)
#'
#' applejuice_uf <- unfold(applejuice) # unfold list into matrix
#' result <- prcomp(applejuice_uf) 
#' # plot PC1 vs PC3 score based on country of production
#' plotScorem(result, ncomp = 4, group = country) 
#'
#' # specify colours
#' plotScorem(result, ncomp = 4, group = country, col = c("black", "grey"))
#' 
#' @importFrom graphics abline legend pairs par plot points text
#' 
#' @export
#' 
plotScorem <-
function(prcompResult, ncomp = 4, group, cex = 1.5, col = NULL, 
         pch = NULL, legendtitle = NULL, ...){
  
  # get information from prcompResult
  score <- prcompResult$x
  N <- length(colnames(score))
  
  # check validity of group
  if (length(group) != dim(score)[1]) {
      stop("The dimension of group and sample do not match. Please check your group variable.")
  }
  if (!is.factor(group) & is.null(attributes(group)$levels)) group <- unclass(as.factor(group))
  
  # number of groups
  numLevels <- nlevels(group)
  
  # color and point type base on group
  col.palette <- generateColor(numLevels, if (!is.null(col))col)
  pch.palette <- generatePoint(numLevels, if (!is.null(pch))pch)
  col <- col.palette[group]
  pch <- pch.palette[group]
  
  # add line for clarity 
  my_line <- function(x,y,...){
    points(x, y, ...)
    abline(v = 0, h = 0, lty = 2, col = "grey39")
  }
  
  # plot
  varlabels <- prcompname(prcompResult, 1:ncomp)
  pairs(score[, 1:ncomp], pch = pch, col = col, 
        lower.panel = my_line, upper.panel = NULL, cex = cex,
        labels = varlabels, ...)
  
  # prepare legend names 
  if (!is.null(attributes(group))) group <- levels(group) # turn unclassed factor back
  
  # add legend
  legend("topright", legend = as.vector(group), pch = pch.palette,
         pt.cex = cex, col =  col.palette, xpd = TRUE, title = legendtitle)
  par(xpd = FALSE)
}
