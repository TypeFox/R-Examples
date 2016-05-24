#' Plot scenario simulation results
#' @description Plots the results like in Study 1 of the Vainik et al. paper an overview of the indicator exclusion results. Starred indicators are excluded in the indicator exclusion procedure. See \code{\link{scenario_sim}} for details. NB! Scenario with no ION violations needs \emph{scenario_plot80()}
#' See  \code{\link{scenario_sim}} for examples

#' @inheritParams scenario_sim
#' @inheritParams scale_sim
#' @param dat simulated data
#' @param multi influences cex of certain plot variables.
#' @param letter assigns plot a letter, useful for combining multiple plots
#' @param ... additional options for axis
#' @param jitter Avoid overalp between lines
#' @encoding utf-8


#' @export
#'
#'
scenario_plot <- function(dat, sizes, n_sim, to_n, tn_n = 8 - to_n, multi = 1, jitter = 0.05, letter = "", ...) {
    old <- options(stringsAsFactors = FALSE)
    # plot n+n simulation results for Study 1
    
    par(mar = c(7, 7, 4, 2))
    
    nlines <- dim(dat)[1]
    
    linenames <- vector()
    
    plot(0, type = "n", axes = F, ann = F, ylim = range(dat), xlim = c(1, length(sizes) + (jitter * nlines)))
    for (i in 1:nlines) {
        ri <- nlines - i + 1
        xpoints <- 1:length(sizes) + (i - 1) * jitter  # set x values for current line
        lines(x = xpoints, y = dat[i, , 1], type = "b", pch = ri, lty = ri, cex = multi, lwd = multi)  # data line
        ciplotter(cix = xpoints, ciy.lo = dat[i, , 2], ciy.hi = dat[i, , 3], eps = 0.04, lty = ri)  # CI-s
        # if(i!=nlines) linenames=c(linenames, paste0(to_n-(i-1),'+',tn_n,' indicators')) else
        # linenames=c(linenames,paste(tn_n, 'indicators')) # use loop to build legend vector
        linenames <- c(linenames, paste0(tn_n, "+", to_n - (i - 1), " indicators"))  # use loop to build legend vector
        
    }
    axis(2, cex.axis = multi, ...)
    axis(1, at = 1:length(sizes), labels = sizes, cex.axis = multi, ...)
    
    # scalename=paste0(letter, 'Scale without ION, ',to_n,'+',tn_n,' indicators, (x',n_sim,')')
    scalename <- paste0(letter, "Scale without ION, ", tn_n, "+", to_n, " indicators")
    
    title(main = scalename, cex.main = multi)
    
    mtext(side = 1, text = "Sample sizes", line = 5, cex = multi * 0.8)
    mtext(side = 2, text = "SONE", line = 5, cex = multi * 0.8)
    
    
    legend(x = "topright", linenames, cex = multi, pch = c(nlines:1), lty = c(nlines:1), lwd = multi, bty = "n")
    on.exit(options(old), add = TRUE)
    
}

#' @describeIn scenario_plot For plotting scenarios where ION is not violated
#' @export
# dat=scen1[[1]] dat=scen1
scenario_plot80 <- function(dat, sizes, n_sim, multi = 1, letter = "", ...) {
    old <- options(stringsAsFactors = FALSE)
    
    # plot 8+0 scale simulation results for Study 1
    par(mar = c(7, 7, 4, 2))
    plot(0, type = "n", axes = F, ann = F, ylim = range(dat), xlim = c(1, length(sizes)))
    lines(dat[, 1], type = "b", pch = 1, lty = 1, cex = multi, lwd = multi)
    axis(2, cex.axis = multi, ...)
    axis(1, at = 1:length(sizes), labels = sizes, cex.axis = multi, ...)
    
    ciplotter(cix = 1:length(sizes), ciy.lo = dat[, 2], ciy.hi = dat[, 3], eps = 0.04)
    
    # scalename=paste0(letter, 'Scale with ION, 8+0 indicators, (x',n_sim,')')
    scalename <- paste0(letter, "Scale with ION, 8+0 indicators")
    
    title(main = scalename, cex.main = multi)
    
    mtext(side = 1, text = "Sample sizes", line = 5, cex = multi * 0.8)
    mtext(side = 2, text = "SONE", line = 5, cex = multi * 0.8)
    
    legend(x = "topright", "8+0 indicators", cex = multi, pch = 1, lty = 1, lwd = multi, bty = "n")
    
    on.exit(options(old), add = TRUE)
    
}


#' Plot confidence intervals
#' @description Used by \code{\link{scenario_plot}}.
#' @param cix values of the x axis
#' @param ciy.lo spread of confidence intervals
#' @param ciy.hi spread of confidence intervals

#' @param eps width of the whiskars
#' @param ... additional base plot arguments
#'

#' @export



ciplotter <- function(cix, ciy.lo, ciy.hi, eps, ...) {
    
    # plot CI-s in base plotter
    segments(cix, ciy.lo, cix, ciy.hi, ...)
    segments(cix - eps, ciy.lo, cix + eps, ciy.lo, ...)
    segments(cix - eps, ciy.hi, cix + eps, ciy.hi, ...)
} 
