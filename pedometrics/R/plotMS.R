#' Model series plot
#' 
#' This function produces a graphical output that allows the examination of the 
#' effect of using different model specifications (design) on the predictive 
#' performance of these models (a model series). It generally is used to access 
#' the results of functions \code{buildMS} and \code{statsMS}, but can be 
#' easily adapted to work with any model structure and performance measure.
#' 
#' @param obj Object of class \code{data.frame}, generally returned by 
#' \code{\link[pedometrics]{statsMS}}, containing a 1) series of performance 
#' statistics of several models, and 2) the design information of each model. 
#' See \sQuote{Details} for more information.
#' @param grid Vector of integer values or character strings indicating the  
#' columns of the \code{data.frame} containing the design data which will be 
#' gridded using the function \code{\link[lattice]{levelplot}}. See 
#' \sQuote{Details} for more information.
#' @param line Character string or integer value indicating which of the 
#' performance statistics (usually calculated by \code{statsMS}) should be 
#' plotted using the function \code{\link[lattice]{xyplot}}. See 
#' \sQuote{Details} for more information.
#' @param ind Integer value indicating for which group of models the mean rank 
#' is to be calculated. See \sQuote{Details} for more information.
#' @param type Vector of character strings indicating some of the effects to be
#' used when plotting the performance statistics using \code{xyplot}. Defaults 
#' to \code{type = c("b", "g")}. See \code{\link[lattice]{panel.xyplot}} 
#' for more information on how to set this argument.
#' @param pch Vector with two integer values specifying the symbols to be used
#' to plot points. The first sets the symbol used to plot the performance
#' statistic, while the second sets the symbol used to plot the mean rank of 
#' the indicator set using argument \code{ind}. Defaults to 
#' \code{pch = c(20, 2)}. See \code{\link[graphics]{points}} for possible values
#' and their interpretation.
#' @param size Numeric value specifying the size of the symbols used for 
#' plotting the mean rank of the indicator set using argument \code{ind}.
#' Defaults to \code{size = 0.5}. See \code{\link[grid]{grid.points}} for more
#' information.
#' @param arrange Character string indicating how the model series should be 
#' arranged, which can be in ascending (\code{asc}) or descending (\code{desc}) 
#' order. Defaults to \code{arrange = "desc"}. See \code{\link[plyr]{arrange}}
#' for more information.
#' @param color Vector defining the colors to be used in the grid produced by 
#' function \code{levelplot}. If \code{NULL}, defaults to 
#' \code{color = cm.colors(n)}, where \code{n} is the number of unique values
#' in the columns defined by argument \code{grid}. See
#' \code{\link[grDevices]{cm.colors}} to see how to use other color palettes.
#' @param xlim Numeric vector of length 2, giving the x coordinates range. If 
#' \code{NULL} (which is the recommended value), defaults to 
#' \code{xlim = c(0.5, dim(obj)[1] + 0.5)}. This is, so far, the optimum range 
#' for adequate plotting.
#' @param ylab Character vector of length 2, giving the y-axis labels. When 
#' \code{obj} is a \code{data.frame} returned by \code{statsMS}, and the
#' performance statistic passed to argument \code{line} is one of those 
#' calculated by \code{statsMS} (\code{"candidates"}, \code{"df"}, 
#' \code{"aic"}, \code{"rmse"}, \code{"nrmse"}, \code{"r2"}, \code{"adj_r2"} or
#' \code{"ADJ_r2"}), the function tries to automatically identify the correct 
#' \code{ylab}.
#' @param xlab Character vector of length 1, giving the x-axis labels. Defaults
#' to \code{xlab = "Model ranking"}.
#' @param at Numeric vector indicating the location of tick marks along the x 
#' axis (in native coordinates).
#' @param ... Other arguments for plotting, although most of these have no been
#' tested. Argument \code{asp}, for example, is not effective since the function
#' automatically identifies the best aspect for plotting based on the dimensions
#' of the design data.
#' 
#' @details
#' This section gives more details about arguments \code{obj}, \code{grid},
#' \code{line}, \code{arrange}, and \code{ind}.
#' \subsection{obj}{
#' The argument \code{obj} usually constitutes a \code{data.frame} returned by
#' \code{statsMS}. However, the user can use any \code{data.frame} object as far 
#' as it contains the two basic units of information needed:
#' \enumerate{
#' \item design data passed with argument \code{grid}
#' \item performance statistic passed with argument \code{line}
#' }
#' }
#' \subsection{grid}{
#' The argument \code{grid} indicates the \emph{design} data which is used to 
#' produce the grid output in the top of the model series plot. By \emph{design} 
#' we mean the data that specify the structure of each model and how they differ 
#' from each other. Suppose that eight linear models were fit using three types 
#' of predictor variables (\code{a}, \code{b}, and \code{c}). Each of these 
#' predictor variables is available in two versions that differ by their 
#' accuracy, where \code{0} means a less accurate predictor variable, while 
#' \code{1} means a more accurate predictor variable. This yields 2^3 = 8 total 
#' possible combinations. The \emph{design} data would be of the following form:
#' \verb{
#' > design
#'   a b c
#' 1 0 0 0
#' 2 0 0 1
#' 3 0 1 0
#' 4 1 0 0
#' 5 0 1 1
#' 6 1 0 1
#' 7 1 1 0
#' 8 1 1 1
#' }
#' }
#' \subsection{line}{
#' The argument \code{line} corresponds to the performance statistic that is
#' used to arrange the models in ascending or descending order, and to produce
#' the line output in the bottom of the model series plot. For example, it can
#' be a series of values of adjusted coefficient of determination, one for each 
#' model:
#' 
#' \verb{
#' adj_r2 <- c(0.87, 0.74, 0.81, 0.85, 0.54, 0.86, 0.90, 0.89)
#' }
#' }
#' \subsection{arrange}{
#' The argument \code{arrange} automatically arranges the model series 
#' according to the performance statistics selected with argument \code{line}.
#' If \code{obj} is a \code{data.frame} returned by \code{statsMS()}, then the
#' function uses standard arranging approaches. For most performance
#' statistics, the models are arranged in descending order. The exception is
#' when \code{"r2"}, \code{"adj_r2"} or \code{"ADJ_r2"} are used, in which case
#' the models are arranged in ascending order. This means that the model with 
#' lowest value appears in the leftmost side of the model series plot, while the
#' models with the highest value appears in the rightmost side of the plot.
#' 
#' \verb{
#' > arrange(obj, adj_r2)
#'   id a b c adj_r2
#' 1  5 1 0 1   0.54
#' 2  2 0 0 1   0.74
#' 3  3 1 0 0   0.81
#' 4  4 0 1 0   0.85
#' 5  6 0 1 1   0.86
#' 6  1 0 0 0   0.87
#' 7  8 1 1 1   0.89
#' 8  7 1 1 0   0.90
#' }
#' 
#' This results suggest that the best performing model is that of \code{id = 7},
#' while the model of \code{id = 5} is the poorest one.
#' }
#' \subsection{ind}{
#' The model series plot allows to see how the design influences model 
#' performance. This is achieved mainly through the use of different colours in
#' the grid output, where each unique value in the \emph{design} data is 
#' represented by a different colour. For the example given above, one could
#' try to see if the models built with the more accurate versions of the
#' predictor variables have a better performance by identifying their relative
#' distribution in the model series plot. The models placed at the 
#' rightmost side of the plot are those with the best performance.
#' 
#' The argument \code{ind} provides another tool to help identifying how the
#' design, more specifically how each variable in the \emph{design} data, 
#' influences model performance. This is done by simply calculating the mean
#' ranking of the models that were built using the updated version of each
#' predictor variable. This very same mean ranking is also used to rank the
#' predictor variables and thus identify which of them is the most important.
#' 
#' After arranging the \code{design} data described above using the adjusted
#' coefficient of determination, the following mean rank is obtained for each
#' predictor variable:
#' 
#' \verb{
#' > rank_center
#'      a    b    c
#' 1 5.75 6.25 5.25
#' }
#' 
#' This result suggests that the best model performance is obtained when using
#' the updated version of the predictor variable \code{b}. In the model series
#' plot, the predictor variable \code{b} appears in the top row, while the 
#' predictor variable \code{c} appears in the bottom row.
#' }
#' @return An object of class \code{"trellis"} consisting of a model series 
#' plot.
#' 
#' @references
#' Deepayan Sarkar (2008). \emph{Lattice: Multivariate Data Visualization with 
#' R.} Springer, New York. ISBN 978-0-387-75968-5.
#' 
#' Roger D. Peng (2008). \emph{A method for visualizing multivariate time series
#' data.} Journal of Statistical Software. v. 25 (Code Snippet), p. 1-17.
#' 
#' Roger D. Peng (2012). \emph{mvtsplot: Multivariate Time Series Plot.} R 
#' package version 1.0-1. \url{http://CRAN.R-project.org/package=mvtsplot}.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
#' @note 
#' Some of the solutions used to build this function were found in the source 
#' code of the R-package \pkg{mvtsplot}. As such, the author of that package, 
#' Roger D. Peng <\email{rpeng@@jhsph.edu}>, is entitled \sQuote{contributors} to
#' the R-package \pkg{pedometrics}.
#' 
#' @section Warning:
#' Use the original functions \code{\link[lattice]{xyplot}} and 
#' \code{\link[lattice]{levelplot}} for higher customization.
#' 
#' @seealso \code{\link[lattice]{levelplot}}, \code{\link[lattice]{xyplot}}, 
#' \code{\link[mvtsplot]{mvtsplot}}.
#' @importFrom stats update
#' @export
#' @examples
#' # This example follows the discussion in section "Details"
#' # Note that the data.frame is created manually
#' id <- c(1:8)
#' design <- data.frame(a = c(0, 0, 1, 0, 1, 0, 1, 1),
#'                      b = c(0, 0, 0, 1, 0, 1, 1, 1),
#'                      c = c(0, 1, 0, 0, 1, 1, 0, 1))
#' adj_r2 <- c(0.87, 0.74, 0.81, 0.85, 0.54, 0.86, 0.90, 0.89)
#' obj <- cbind(id, design, adj_r2)
#' p <- plotMS(obj, grid = c(2:4), line = "adj_r2", ind = 1, 
#'             color = c("lightyellow", "palegreen"),
#'             main = "Model Series Plot")
#' print(p)
#' 
#' @keywords hplot
#' 
# FUNCTION #####################################################################
plotMS <-
  function (obj, grid, line, ind, type = c("b", "g"), pch = c(20, 2),
            size = 0.5, arrange = "desc", color = NULL, 
            xlim = NULL, ylab = NULL, xlab = NULL, at = NULL, ...) {
    
    # Check if suggested packages are installed
    pkg <- c("grDevices", "grid", "plyr")
    id <- !sapply(pkg, requireNamespace, quietly = TRUE)
    if (any(id)) {
      pkg <- paste(pkg[which(id)], collapse = " ")
      stop(paste("Package(s) needed for this function to work but not",
                 "installed: ", pkg, sep = ""), call. = FALSE)
    }
    
    if (missing(obj)) {
      stop("<obj> is a mandatory argument")
    }
    if (missing(grid)) {
      stop("<grid> is a mandatory argument")
    }
    if (missing(line)) {
      stop("<line> is a mandatory argument")
    }
    if (missing(ind)) {
      stop("<ind> is a mandatory argument")
    }
    if (class(obj) != "data.frame") {
      stop("<obj> should be of class data.frame")
    }
    if (!any(class(grid) == c("integer", "character", "numeric"))) {
      stop("<grid> should be an integer value or a character string")
    }
    if (!any(class(line) == c("integer", "character", "numeric"))) {
      stop("<line> should be an integer value or a character string")
    }
    if (!any(class(ind) == c("integer", "numeric")) || round(ind) != ind) {
      stop("<ind> should be an integer value")
    }
    if (any(class(line) == c("integer", "numeric"))) {
      nam0 <- c("candidates", "df", "aic", "rmse", "nrmse", "r2", "adj_r2", 
                "ADJ_r2")
      nam1 <- colnames(obj)[line]
      if (!any(colnames(obj)[line] == nam0)) {
        stop(paste("<ylab> should be provided for performance statistics <",
                   nam1, ">",  sep = ""))
      }
    }
    if (!missing(xlab)) {
      if (length(xlab) != 1) {
        stop("<xlab> should have length equal to 1")
      }
    }
    if (!missing(ylab)) {
      if (length(ylab) != 2) {
        stop("<ylab> should have length equal to 2")
      }
    }
    if (length(type) != 2) {
      stop("<type> should have length equal to 2")
    }
    if (length(pch) != 2) {
      stop("<pch> should have length equal to 2")
    }
    # prepare data #############################################################
    if (class(line) == "numeric") {
      line <- colnames(obj)[line]
    }
    if (any(line == c("r2", "adj_r2", "ADJ_r2"))) {
      obj <- plyr::arrange(obj, plyr::desc(obj[, line]))
    } else {
      obj <- plyr::arrange(obj, obj[, line])
    }
    grid <- as.matrix(obj[, grid])
    x <- seq(1, dim(obj)[1], 1)
    y <- as.numeric(obj[, line])
    if (missing(at)) {
      if (max(x) < 100) {
        m <- round(max(x) / 10) * 10
        at <- c(1, seq(5, m, 5))
      } else {
        m <- round(max(x) / 10) * 10
        at <- c(1, seq(10, m, by = 10))
      }
    }
    if (missing(color)) {
      color <- grDevices::cm.colors(length(unique(as.numeric(grid))))
    }
    if (missing(xlim)) {
      xlim <- c(0.5, dim(obj)[1] + 0.5)
    }
    if (missing(xlab)) {
      xlab <- "Model ranking"
    }
    if (missing(ylab)){
      if (class(line) == "numeric") {
        line <- colnames(obj)[line]
      }
      if (line == "candidates") {
        yl <- "Candidate predictors"
      }
      if (line == "df") {
        yl <- "Degrees of freedom"
      }
      if (line == "aic") {
        yl <- "AIC"
      }
      if (line == "rmse") {
        yl <- "RMSE"
      }
      if (line == "nrmse") {
        yl <- "NRMSE"
      }
      if (line == "r2") {
        yl <- expression(paste(R^2, sep = ''))
      }
      if (any(line == c("adj_r2", "ADJ_r2"))) {
        yl <- expression(paste('Adjusted ',R^2, sep = ''))
      }
      ylab <- list(c(yl, "Design"))
    }
    rank_center <- rep(NA, dim(grid)[2])
    for (i in 1:length(rank_center)) {
      rank_center[i] <- 
        mean(cbind(x, grid)[, 1][which(cbind(x, grid)[, i + 1] == ind)])
    }
    grid <- grid[, order(rank_center, decreasing = TRUE)]
    p1 <- lattice::xyplot(
      y ~ x, xlim = rev(grDevices::extendrange(xlim, f = 0)), type = type, 
      pch = pch[1], scales = list(y = list(rot = 0), x = list(at = at)))
    p2 <- lattice::levelplot(
      grid, colorkey = FALSE, xlim = rev(grDevices::extendrange(xlim, f = 0)),
      col.regions = color, scales = list(y = list(rot = 90)),
      panel = function (...) {
        lattice::panel.levelplot(...)
        grid::grid.points(x = sort(rank_center, decreasing = TRUE), 
                    seq(1, dim(grid)[2], 1),
                    pch = pch[2], size = grid::unit(size, "char"))
        })
    
    # Print plot
    update(c(p1, p2), layout = c(1, 2), xlab = xlab, 
           ylab = ylab, aspect = c((dim(grid)[2] * 2) / dim(grid)[1]),
           par.settings = list(layout.heights = list(panel = c(0.5, 0.5))), ...)
    }

