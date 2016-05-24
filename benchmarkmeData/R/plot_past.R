nice_palette = function(){
  alpha =150
  palette(c(rgb(85,130,169, alpha=alpha, maxColorValue=255),
            rgb(200,79,178, alpha=alpha,maxColorValue=255), 
            rgb(105,147,45, alpha=alpha, maxColorValue=255),
            rgb(204,74,83, alpha=alpha, maxColorValue=255),
            rgb(183,110,39, alpha=alpha, maxColorValue=255),
            rgb(131,108,192, alpha=alpha, maxColorValue=255)))
}


#' Plot past results
#' 
#' Plot the previous benchmarks. This function creates two figures.
#' \itemize{
#' \item Figure 1: Total benchmark time over all benchmarks (in seconds) on the y-axis.
#' \item Figure 2: Relative time (compared to the smallest benchmark).
#' }
#' The data set used is \code{data(past_results)}.
#' @param test_group One of "prog", "matrix_fun", "matrix_cal", "read", "write". Default it \code{prog}.
#' @param byte_optimize Default \code{NULL}. The default behaviour is to plot all results.
#' To plot only the byte optimized results, set to \code{TRUE}, otherwise \code{FALSE}.
#' @param blas_optimize Default \code{NULL}. The default behaviour is to plot all results.
#' To plot only the BLASS optimized results, set to \code{TRUE}, otherwise \code{FALSE}.
#' @param log By default the y axis is plotted on the log scale. To change, set the 
#' the argument equal to the empty parameter string, \code{""}.
#' @examples 
#' ## Plot non byte optimize code
#' plot_past(byte_optimize=FALSE)
#' @importFrom graphics abline grid par plot points legend
#' @importFrom grDevices palette rgb
#' @importFrom utils data
#' @importFrom stats aggregate
#' @export
#' @examples 
#' ## Plot all past results
#' plot_past()
#' 
#' ## Plot the blas_optimized results
#' plot_past(blas_optimize=TRUE)
plot_past = function(test_group=c("prog", "matrix_fun", "matrix_cal", 
                                  "read", "write"), 
                     byte_optimize=NULL, blas_optimize=NULL,
                     log="y") {
  

  results = select_results(test_group, byte_optimize, blas_optimize)
  
  
  ## Arrange plot colours and layout
  op = par(mar=c(3,3,2,1), 
           mgp=c(2,0.4,0), tck=-.01,
           cex.axis=0.8, las=1, mfrow=c(1,2)) 
  old_pal = palette()
  on.exit({palette(old_pal); par(op)})
  nice_palette()
  
  
  ymin = min(results$time)
  ymax = max(results$time)
  plot(results$time, xlab="Rank", ylab="Total timing (secs)", 
       ylim=c(ymin, ymax), xlim=c(1, nrow(results)+1), cex=0.9,
       panel.first=grid(), log=log, pch=21, bg=as.numeric(results$test_group))
  
  if(any(c("read5", "write5") %in% results$test_group)) {
    legend("topleft", pt.bg=1:4, 
           pch=21, legend = levels(results$test_group), 
           bty="n")
  }
  
  ## Relative timings  
  fastest = min(results$time)
  ymax = ymax/fastest
  plot(results$time/fastest, xlab="Rank", ylab="Relative timing", 
       ylim=c(1, ymax), xlim=c(1, nrow(results)+1), cex=0.9,
       panel.first=grid(), log=log, pch=21, bg=as.numeric(results$test_group))
  abline(h=1, lty=3)
  
  invisible(results)
}
