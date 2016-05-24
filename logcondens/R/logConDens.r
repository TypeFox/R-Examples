logConDens <- function(x, xgrid = NULL, smoothed = TRUE, print = FALSE, gam = NULL, xs = NULL){

    res1 <- activeSetLogCon(x, xgrid = xgrid, print = print)
    if (identical(smoothed, FALSE)){res <- c(res1)}

    if (identical(smoothed, TRUE)){
    
        if (identical(xs, NULL)){
            r <- diff(range(x))
            xs <- seq(min(x) - 0.1 * r, max(x) + 0.1 * r, length = 500)
        }
    
        smo <- evaluateLogConDens(xs, res1, which = 4:5, gam = gam, print = print)
        f.smoothed <- smo[, "smooth.density"]
        F.smoothed <- smo[, "smooth.CDF"]
        mode <- xs[f.smoothed == max(f.smoothed)]
        
        res2 <- list("f.smoothed" = f.smoothed, "F.smoothed" = F.smoothed, "gam" = gam, "xs" = xs, "mode" = mode)    
        res <- c(res1, res2)
        }

    res$smoothed <- smoothed

    class(res) <- "dlc"
    return(res)
}
