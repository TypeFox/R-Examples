#' Plot Panel
#' 
#' The \code{plot_panel} function takes objects of the class
#' \code{\linkS4class{adpcr}} to enable customizable graphical representations
#' of a chamber-based digital PCR experiments (e.g., Digital Array (R) IFCs
#' (integrated fluidic circuits) of the BioMark (R) and EP1 (R)).
#' 
#' @details Currently, only objects containing \code{tnp} data can be plotted as
#' a whole. For the any other type of the \code{adpcr} data, only just one column 
#' of data (one panel) can be plotted at the same time (see Examples how easily 
#' plot multipanel objects). Moreover the object must contain fluorescence 
#' intensities or exact number of molecules or
#' the positive hits derived from the Cq values for each well. The Cq values
#' can be obtained by custom made functions (see example in
#' \code{\link{dpcr_density}})) or the yet to implement "qpcr_analyser function
#' from the dpcR package.
#' 
#' If the \code{col} argument has length one, a color is assigned for each
#' interval of the input, with the brightest colors for the lowest values.
#' 
#' @inheritParams adpcr2panel
#' @param col A single color or vector of colors for each level of input.
#' @param legend If \code{TRUE}, a built-in legend is added to the plot.
#' @param half If \code{left} or \code{right}, every well is represented only
#' by the adequate half of the rectangle.
#' @param plot \code{"logical"}, if \code{FALSE}, only plot data is returned
#' invisibly.
#' @param ... Arguments to be passed to \code{plot} function.
#' @return Invisibly returns two sets of coordinates of each microfluidic well
#' as per \code{\link{calc_coordinates}}:
#' \code{coords} is a list of coordinates suitable for usage with functions from
#' \code{\link{graphics}} package. The second element is a data frame of coordinates 
#' useful for users utilizing ggplot2 package.
#' @author Michal Burdukiewicz, Stefan Roediger.
#' @seealso \code{\link{extract_dpcr}} - extract experiments.
#' \code{\link{adpcr2panel}} - convert \code{\linkS4class{adpcr}} object to arrays.
#' @keywords hplot
#' @examples
#' 
#' # Create a sample dPCR experiment with 765 elements (~> virtual compartments)   
#' # of target molecule copies per compartment as integer numbers (0,1,2)
#' ttest <- sim_adpcr(m = 400, n = 765, times = 20, pos_sums = FALSE, 
#'                    n_panels = 1)
#' # Plot the dPCR experiment results with default settings
#' plot_panel(ttest)
#' 
#' #do it without breaks
#' plot_panel(ttest, use_breaks = FALSE)
#' 
#' # Apply a binary color code with blue as positive
#' slot(ttest, "breaks") <- c(0, 2, 4)
#' plot_panel(ttest, col = "blue")
#' 
#' # Apply a two color code for number of copies per compartment
#' plot_panel(ttest, col = c("blue", "red"))
#' 
#' 
#' 
#' # supply customized breaks and compare
#' par(mfcol = c(2, 1))
#' plot_panel(ttest)
#' slot(ttest, "breaks") <- c(0, 1, 2, (max(slot(ttest, "breaks")) + 1))
#' plot_panel(ttest)
#' par(mfcol = c(1, 1))
#' 
#' # plot few panels
#' ttest2 <- sim_adpcr(m = 400, n = 765, times = 40, pos_sums = FALSE, 
#'                     n_panels = 4)
#' par(mfcol = c(2, 2))
#' four_panels <- lapply(1:ncol(ttest2), function(i) 
#'        plot_panel(extract_dpcr(ttest2, i), legend = FALSE, 
#'          main = paste("Panel", LETTERS[i], sep = " ")))
#' par(mfcol = c(1, 1))
#' 
#' # two different channels 
#' plot_panel(extract_dpcr(ttest2, 1), legend = FALSE, 
#'            half = "left")
#' par(new = TRUE)
#' plot_panel(extract_dpcr(ttest2, 2), col = "blue", 
#'            legend = FALSE, half = "right")
#' 
#' # plot two panels with every well as only the half of the rectangle
#' ttest3 <- sim_adpcr(m = 400, n = 765, times = 40, pos_sums = FALSE, 
#'                     n_panels = 2)
#' par(mfcol = c(1, 2))
#' two_panels <- lapply(1:ncol(ttest3), function(i) 
#'        plot_panel(extract_dpcr(ttest3, i), legend = FALSE, 
#'          main = paste("Panel", LETTERS[i], sep = " ")))
#' par(mfcol = c(1, 1))
#' 
#' @export plot_panel
plot_panel <- function(input, use_breaks = TRUE, col = "red", legend = TRUE, 
                       half = "none", plot = TRUE, ...) {  
  
  array <- adpcr2panel(input, use_breaks = use_breaks)
  
  if(length(array) > 1) 
    warning("Only the first array will be processed.")

  array <- array[[1]]
  all_coords <- calc_coordinates(array, half = half)

  if(plot) {
    nx_a <- ncol(array) 
    ny_a <- nrow(array)
    
    cutted_input <- factor(array)
    cols <- cutted_input
    ncols <- nlevels(cutted_input)
    if (length(col) == 1) {   
      levels(cols) <- sapply(0:ncols/ncols, function(x) 
        adjustcolor(col, alpha.f = x))
    } else {
      if (length(col) != ncols) {
        stop("The vector of colors must have length equal to the number of levels of 
             the input.")    
      }
      levels(cols) <- col
    }
    
    plot(NA, NA, xlim = c(1, nx_a), ylim = c(1, ny_a), axes = FALSE, xlab = "", 
         ylab = "", ...)
    
    if (legend)
      legend(x = -0.085 * nx_a, 
             y = ny_a/1.6, 
             legend = levels(cutted_input),
             fill = levels(cols), 
             bty = "n", 
             xpd = TRUE, 
             x.intersp = 0.5)
    
    cols <- as.character(cols)
    args <- lapply(1L:length(all_coords[["coords"]]), function(i) 
      c(all_coords[["coords"]][[i]], list(col = cols[i])))
    
    sapply(1L:length(input), function(i) 
      do.call(rect, args[[i]]))
  }
  
  invisible(list(coords = all_coords[["coords"]], ggplot_coords = all_coords[["ggplot_coords"]]))
}



