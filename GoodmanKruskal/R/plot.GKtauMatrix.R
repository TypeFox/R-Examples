#' Plot method for DataRobot S3 objects of class GKtauMatrix
#'
#' Method for R's generic plot function for DataRobot S3 objects
#' of class GKtauMatrix.  This function generates an array of
#' Goodman-Kruskal tau association measures as described under
#' Details.  Note that, in general, this matrix is asymmetric.
#'
#' This function calls the corrplot function from the corrplot
#' package to generate an array of correlation plots from the matrix
#' of Goodman-Kruskal tau measures returned by the GKtauDataframe
#' function.  The off-diagonal elements of this array contain
#' ellipses whose shape is determined by the corresponding tau
#' value, and the diagonal elements display the number of unique
#' levels for each variable in the plot.
#'
#' This plot may be rendered either in color (the default, obtained
#' by specifying colorPlot as TRUE) or black-and-white.  In color
#' plots, the color of the text for the correlation values is set
#' by the corrColors parameter.  The default value for this parameter
#' is NULL, which causes the function to use the color vector
#' rainbow(n) where n is the number of rows and columns in the
#' GKtauMatrix object x.  The background color used to fill in
#' each ellipse is specified by bhe backgroundColor parameter,
#' and the text for the diagonal entries is determined by the
#' diagColor parameter.  In cases where the default choices make
#' the correlation values difficult to read, a useful alternative
#' is to specify corrColors = "blue".
#'
#' @param x S3 object of class GKtauMatrix to be plotted.
#' @param y Not used; included for conformance with plot() generic
#' function parameter requirements.
#' @param colorPlot Logical variable indicating whether to
#' generate a color plot (the default, for colorPlot = TRUE)
#' or a black-and-white plot.
#' @param corrColors Character vector giving the color names
#' for the correlation values printed on the plot; default
#' value is NULL, causing rainbow(n) to be used, where n is
#' the number of rows and columns in the matrix x.
#' @param backgroundColor Character variable naming the background
#' color used for the correlation ellipses in the plot.
#' @param diagColor Character variable naming the color of the
#' text used to display the number of levels per variable along
#' the diagonal of the correlation matrix plot.
#' @param diagSize Numeric scale factor to adjust the text size for
#' the number of levels displayed on the diagonal of the plot array.
#' @param \dots Not used; included for conformance with plot() generic
#' function parameter requirements.
#' @return None.  This function is called for its side-effect of
#' generating a plot.
#' @author Ron Pearson
#' @export
#'
plot.GKtauMatrix <- function(x, y, colorPlot = TRUE,
                             corrColors = NULL,
                             backgroundColor = "gray",
                             diagColor = "black",
                             diagSize = 1, ...){
  #
  GKclip <- pmin(x, 1)
  n <- ncol(x)
  #
  if (colorPlot){
    if (is.null(corrColors)){
      corrColors <- rainbow(n)
    }
    corrplot::corrplot(GKclip, method = "ellipse", outline = TRUE,
             col = backgroundColor, diag = FALSE, cl.pos = "n",
             addCoef.col = corrColors, tl.srt = 45,
             tl.col = "black")
  } else {
    corrplot::corrplot(GKclip, method = "ellipse", outline = TRUE,
             col = "transparent", diag = FALSE, cl.pos = "n",
             addCoef.col = "black", tl.srt = 45,
             tl.col = "black")
  }
  for (i in 1:n){
    topLine <- paste("K =", x[i,i])
    text(i, n + 1 - i, topLine, col = diagColor, cex = diagSize)
  }
}
