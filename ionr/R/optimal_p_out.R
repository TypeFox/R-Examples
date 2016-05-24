#' Table and plot the SONE values
#'
#' @description Takes max and min scenarious and produces a table and optionally a plot.
#'  See \code{\link{scenario_sim}} or \code{\link{optimal_p}}.
#' @inheritParams optimal_p
#' @param scenario_max SONE data from \code{\link{scenario_sim}} output
#' @param scenario_min SONE data from \code{\link{scenario_sim}} output
#' @param multi influences cex of certain plot variables. Defaults to 1
#' @encoding utf-8

#' @examples
#' set.seed(466)
#' sizes=c(500,1000)
#' n_sim=50  #  make bigger for more accurate estimates..
#' to_n=8
#' cor_to_outcome=0.25
#' ptm <- proc.time()  # timing
#' # takes a few seconds..
#' scen1=scenario_sim(sizes=sizes,n_sim=n_sim,to_n=to_n, cor_to_outcome=cor_to_outcome)
#' proc.time() - ptm
#' ptm <- proc.time()
#' # A scenario with 3 out of 8 items relating to outcome, 3 different samples
#' to_n=3
#' scen2=scenario_sim(sizes=sizes,n_sim=n_sim,to_n=to_n, cor_to_outcome=cor_to_outcome)
#' proc.time() - ptm
#'
#' optimal_p_out(scen1[[1]],scen2[[1]],sizes = sizes,n_sim=n_sim,to_min = to_n, plot='yes', multi=1)
#'
#' # Should be equivalent. Some variation can be expected when n_sim is below 1000
#' ptm <- proc.time()
#' a=optimal_p(sizes=sizes, n_sim=n_sim, n_indicators=8, plotting='yes', cor_to_outcome=cor_to_outcome)
#' proc.time() - ptm
#' print(a[[1]])
#' @export

optimal_p_out <- function(scenario_max, scenario_min, sizes, n_sim, to_min, plotting = "", multi = 1) {
    
    # single sizes values need to be handled differently
    if (length(sizes) == 1) {
        if (plotting != "") {
            plotting = ""
            print("Plotting disabled, as 'sizes' is a single value")
        }
        
        crithi <- scenario_max[1]
        # p value below which 1+n_indicators - to_min indicator is not detected
        critlo <- scenario_min[to_min, 1]
        
    } else {
        crithi <- scenario_max[, 1]
        # p value below which 1+n_indicators - to_min indicator is not detected
        critlo <- scenario_min[to_min, , 1]
    }
    # see that max > min
    critdelta <- crithi - critlo
    
    # get the optimal values
    crit <- crithi  # temp variable
    
    for (i in 1:length(crithi)) {
        crit[i] <- psych::geometric.mean(c(crithi[i], critlo[i]), na.rm = TRUE)
    }
    
    
    
    # if highest pvalue is lower than lowest p value, then we cannot really find an optimal p value
    crit[critdelta < 0] <- NA
    
    out <- matrix(ncol = length(sizes), nrow = 3)
    dimnames(out) <- list(c("optimal", "max", "min"), sizes)
    out[1, ] <- crit
    out[2, ] <- crithi
    out[3, ] <- critlo
    
    if (plotting == "yes" | plotting == "file") {
        old <- options(stringsAsFactors = FALSE)
        
        # plot the p values multi=1.5
        if (plotting == "file") 
            tiff(paste0("optimal p_x", n_sim, "_", Sys.Date(), ".tiff"), res = 300, width = 30, height = 20, units = "cm", 
                compression = "lzw", pointsize = 15)
        # par(mfrow = c(1, 1), mar = c(5, 6, 4, 2))
        
        
        plot(0, type = "n", axes = F, ann = F, ylim = range(critlo, crithi), xlim = c(1, length(sizes)))
        lines(critlo, type = "b", pch = 3, lty = 3, cex = multi, lwd = multi)
        lines(crithi, type = "b", pch = 4, lty = 2, cex = multi, lwd = multi)
        lines(x = 1:length(sizes), y = crit, type = "b", pch = 1, lty = 1, cex = multi, lwd = multi)  # ignore crit[1],as critdelta is negative
        axis(2, cex.axis = multi)
        axis(1, at = 1:length(sizes), labels = sizes, cex.axis = multi)
        
        # linenames <- c('Minimum criteria', 'Maximum criteria', 'Optimal criteria') legend(x = 'topright', linenames, cex
        # = multi, pch = c(3, 4, 1), lty = c(3:1), lwd = multi, bty = 'n') Match the sequence of the lines in plot
        linenames <- c("Maximum criteria", "Optimal criteria", "Minimum criteria")
        legend(x = "topright", linenames, cex = multi, pch = c(4, 1, 3), lty = c(2, 1, 3), lwd = multi, bty = "n")
        
        
        # title(main='Optimal p criteria for sample sizes', xlab='sample sizes',ylab='Minimum 'significance of indicator
        # exclusion'', cex.main=multi, cex.lab=multi)
        title(xlab = "Sample sizes", ylab = "SONE", cex.main = multi, cex.lab = multi)
        
        if (plotting == "file") 
            dev.off()
        on.exit(options(old), add = TRUE)
        
    }
    
    return(out)
    
    
}
#  
