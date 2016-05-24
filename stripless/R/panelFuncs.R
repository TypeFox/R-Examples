#' strucplot Panel Functions
#' 
#' Panel functions for strucplot.
#' 
#' @details panel.bars A wrapper for \code{\link[lattice]{panel.barchart}} that 
#'   plots bars that summarize the responses at each setting of the 
#'   conditioning variables. By default, bars are vertical, but setting the
#'   optional parameter, \code{horizontal} to \code{TRUE} plots horizontal bars.
#'   
#' @param ... Arguments to barchart. See \code{\link[lattice]{panel.barchart}}.
#' 
#' @param summaryFUN The function that summarizes replicated response values.
#' Default = fuction(x) mean(x, na.rm = TRUE).
#' @param col Color of the bars. Default = "darkblue".
#' @param grid Should an appropriate background grid be plotted? Default = TRUE.
#' @param col.grid The background grid color. Default = "lightgray".
#' @examples
#' # A half fraction of a 2^5 full factorial with pseudo-replicate responses
#' # at each design point,
#' 
#' # Build the design matrix
#' x <- c(-1,1)
#' ff <- expand.grid(x,x,x,x)
#' ff[[5]] <- do.call(mapply,c(FUN=prod,ff))
#' ff <- ff[rep(1:16,e=2),] ## replicates each row twice
#' names(ff) <- LETTERS[1:5]
#' 
#' # Add a column for the response
#' ff$y <-c(155.5, 154.8, 158.4, 156.2, 154.8, 152.4, 159.7, 155.5, 161.8, 
#' 159.7, 159, 158.4, 159.7, 157.7, 161.8, 158, 155.9, 151.7, 159, 
#' 158, 154.1, 156.9, 158.4, 158.4, 159, 154.8, 158.4, 156.2, 161.1, 
#' 156.9, 162.6, 159)
#' 
#' # Plot using panel.bars
#' strucplot(~ y|., data = ff, panel = panel.bars)
#' 
#' # It is often useful to plot the bars the other way, too
#' strucplot(~ y|., data = ff, panel = panel.bars, horizontal = TRUE)
#' 

panel.bars <- function(..., summaryFUN = function(x)mean(x,na.rm=TRUE),
                       col = "darkblue",grid=TRUE, col.grid="lightgray")
{
 dots <- list(...)
 hz <- dots[["horizontal"]]
 if(is.null(hz) || !is.logical(hz)) {hz <- dots[["horizontal"]] <-FALSE}
 if(grid) do.call(panel.grid,c(list(h = -!hz, v = -hz, col.line = col.grid), dots))
 x <- dots[["x"]]; y <- dots[["y"]]
 if(any(vapply(list(x,y),function(z)!length(z), TRUE))
     ||anyNA(c(x,y)))  return()   ## empty plot
 else {
   dots[["col"]] <- col
   if(hz) dots[["x"]] <- ave(x,y, FUN=summaryFUN)
   else dots[["y"]] <- ave(y, x, FUN=summaryFUN)
  }
 do.call(panel.barchart,dots)
}
 

