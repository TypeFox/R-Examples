#' @title Run Evanno Method on STRUCTURE Results
#' @description Calculate first and second order rates of changes of LnPr(K) 
#'   from STRUCTURE results based on Evanno et al. 2005.
#' 
#' @param sr output from a call to \code{\link{structure}}.
#' @param plot logical. Generate a plot of Evanno metrics.
#' 
#' @return a data.frame with Evanno log-likelihood metrics for each value of K.
#' 
#' @references Evanno, G., Regnaut, S., and J. Goudet. 2005. Detecting the 
#'   number of clusters of individuals using the software STRUCTURE: a 
#'   simulation study. Molecular Ecology 14:2611-2620.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{structure}}
#' 
#' @importFrom stats sd
#' @importFrom graphics par layout segments axis box text
#' @export
 
evanno <- function(sr, plot = TRUE) {
  if(!"structure.result" %in% class(sr)) {
    stop("'sr' is not a result from 'structure.run'.")
  }
  k.tbl <- table(sapply(sr, function(x) x$summary["k"]))
  if(length(k.tbl) < 3) stop("must have at least two values of k.")
  
  # collect summary statistics
  sr.smry <- t(sapply(sr, function(x) x$summary))
  
  # calculate mean and sd of LnPr(K)
  ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[, "k"], mean)
  sd.ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[, "k"], sd)
  
  # Ln'(K)
  ln.pk <- diff(ln.k)
  # Ln''(K)
  ln.ppk <- abs(diff(ln.pk))
  # Delta-K
  delta.k <- sapply(2:(length(ln.k) - 1), function(i) {
    abs(ln.k[i + 1] - (2 * ln.k[i]) + ln.k[i - 1]) / sd.ln.k[i]
  })
  
  result <- data.frame(
    k = as.numeric(names(ln.k)),
    reps = as.numeric(table(sr.smry[, "k"])),
    mean.ln.k = as.numeric(ln.k),
    sd.ln.k = as.numeric(sd.ln.k),
    ln.pk = c(NA, ln.pk),
    ln.ppk = c(NA, ln.ppk, NA),
    delta.k = c(NA, delta.k,  NA)
  )
  rownames(result) <- NULL
  
  if(plot) {
    xlim <- range(c(0, result$k))
    tick.spacing <- if(max(result$k) < 10) {
      1
    } else if(max(result$k) < 20) {
      2
    } else 5
    xticks <- seq(0, max(result$k), tick.spacing)
    
    op <- par(mar = c(4, 4, 1, 1) + 0.1)
    layout(matrix(1:4, nrow = 2, byrow = TRUE))
    
    plot.func <- function(y, ylab, sd = NULL) {
      ylim <- if(is.null(sd)) {
        range(y, na.rm = TRUE) 
      } else {
        sd[is.na(sd)] <- 0
        range(c(y, y + sd, y - sd), na.rm = TRUE)
      }     
      plot(result$k, y, xlim = xlim, ylim = ylim, type = "b", 
           xlab = "K", ylab = ylab, axes = F, 
           bty = "l", pch = 19, bg = "black")
      if(!is.null(sd)) segments(result$k, y + sd, result$k, y - sd)
      axis(1, at = xticks[-1])
      axis(2)
      box(bty = "l")
    }
    
    plot.func(result$mean.ln.k, "mean LnP(K)", sd = result$sd.ln.k)
    plot.func(result$ln.pk, "LnP'(K)")
    plot.func(result$ln.ppk, "LnP''(K)")
    if(!all(is.na(result$delta.k))) {
      plot.func(result$delta.k, "Delta(K)")
    } else {
      plot(0, type = "n", ann = FALSE, axes = FALSE)
      text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), "N/A")
    }
    
    layout(matrix(1))
    par(op)
  } 
  
  result
}