

gplot <- function(x, h, xlim = range(object$x, na.rm = TRUE), cex = 0.7, nplot = 3, type = c("base", "ggplot"), ...) {
    object <- x
    cutpoint <- getCutpoint(object)
    
    ## bandwidth: use Ruppert, Sheather and Wand (KernSmooth:::dpill)
    if (missing(h)) {
        if (!all(xlim == range(object$x, na.rm = TRUE))) {
            object <- subset(object, object$x > min(xlim) & object$x < max(xlim))
        }
        h <- rdd_bw_rsw(object)
        if (is_even(nplot)) {
            se <- seq(from = 1 - (sum(1:nplot < (nplot/2))) * 0.2, to = 1 + (sum(1:nplot > (nplot/2))) * 0.2, by = 0.2)
        } else {
            se <- seq(from = 1 - floor(nplot/2) * 0.2, to = 1 + floor(nplot/2) * 0.2, by = 0.2)
        }
        hs <- if (nplot == 1) 
            h else se * h
    } else {
        if (length(h) == 1) {
            if (is_even(nplot)) {
                se <- seq(from = 1 - (sum(1:nplot < (nplot/2))) * 0.2, to = 1 + (sum(1:nplot > (nplot/2))) * 0.2, by = 0.2)
            } else {
                se <- seq(from = 1 - floor(nplot/2) * 0.2, to = 1 + floor(nplot/2) * 0.2, by = 0.2)
            }
            hs <- if (nplot == 1) 
                h else se * h
        } else {
            if (length(h == nplot)) {
                hs <- h
            } else {
                stop("Length of h should be either one or equal to nplot (", nplot, ")")
            }
        }
    }
    
    
    
    
    ## plot
    if (type == "base") {
        par_orig <- par()
        par(mfrow = c(nplot, 1))
        for (i in 1:nplot) {
            plotBin(x = object$x, y = object$y, cutpoint = cutpoint, h = hs[i], xlim = xlim, cex = cex)
        }
    } else {
        
        plotBin_out <- plotBin(x = object$x, y = object$y, cutpoint = cutpoint, h = hs[1], xlim = xlim, cex = cex, plot = FALSE)
        plotBin_out$h <- rep(hs[1], nrow(plotBin_out))
        for (i in 2:nplot) {
            new <- plotBin(x = object$x, y = object$y, cutpoint = cutpoint, h = hs[i], xlim = xlim, cex = cex)
            new$h <- rep(hs[i], nrow(new))
            plotBin_out <- rbind(plotBin_out, new)
        }
        
        plotBin_out$h <- round(plotBin_out$h, 4)
        qplot(x = x, y = y, data = plotBin_out) + facet_grid(h ~ .)
        
    }
    
} 
