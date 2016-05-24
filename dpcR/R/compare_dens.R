#' Plot and Compare Densities
#' 
#' Plot empirical and theoretical density of the result of the digital PCR
#' experiment.
#' 
#' 
#' @param input object of class \code{\linkS4class{dpcr}} containing only one run.
#' @param moments logical, if \code{TRUE}, both theoretical and empirical
#' moments are printed on the plot.
#' @param \dots other arguments passed to the \code{plot} function.
#' @author Michal Burdukiewcz.
#' @seealso \link{moments} is used to calculate moments of Poisson distribution.
#' @keywords density digital PCR empirical
#' @examples
#' 
#' adpcr_big <- sim_adpcr(m = 35, n = 40, times = 50, pos_sums = FALSE, n_panels = 1)
#' compare_dens(adpcr_big, moments = TRUE)
#' 
#' @export compare_dens
compare_dens <- function(input, moments = TRUE, ...) {  
  #moments() checks class and so on
  
  if (ncol(input) > 1)
    stop("Input must contain only one experiment.")    
  
  all_moms <- moments(input)
  lambda <- all_moms[all_moms[["moment"]] == "mean", "value"][1]
  xup <- max(input)
  data <- table(factor(input, levels = 0L:xup))
  bars <- calc_bars(data)
  theor <- dpois(0L:xup, lambda)*slot(input, "n")
  ytop <- ifelse(max(theor) >= max(data), max(theor), max(data))
  
  plot(NA, NA, xlim = c(-0.5, xup + 0.5), ylim = c(-0, ytop), 
       xlab = "Number of molecules", ylab = "Counts", ...)
  
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = adjustcolor("grey", alpha.f = 0.30))
  axis(2, tck = 1, col.ticks = "white", labels = FALSE)
  
  apply(bars, 1, function(x) 
    rect(x[1], x[2], x[3], x[4]))
  #   axis(4, at = theor, labels = 0L:xup, tck = 1, lty = "dotted", 
  #        col.ticks = "darkgrey")
  #   mtext("Theoretical counts", side = 4, line = 2) 
  sapply(0L:xup, function(x) 
    lines(c(x, x), c(0, theor[x + 1]), lty = "dotted", col = "grey12", lwd = 2))
  
  if (moments) {
    labels <- as.character(unique(all_moms[["moment"]]))
    sapply(1L:4, function(i) {
      text(0.83*xup, (98 - 5*i)/100*ytop, paste0(labels[i], ":"), pos = 2)
      text(0.87*xup, (98 - 5*i)/100*ytop, round(all_moms[all_moms[["method"]] == "theoretical", "value"][i], 4))
      text(0.99*xup, (98 - 5*i)/100*ytop, round(all_moms[all_moms[["method"]] == "empirical", "value"][i], 4))
    })
    text(0.87*xup, 0.99*ytop, "Theoretical", pos = 1)
    text(0.99*xup, 0.99*ytop, "Empirical", pos = 1)
  }
}
