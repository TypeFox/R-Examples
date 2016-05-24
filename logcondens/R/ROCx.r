ROCx <- function(x, res0, res1, smooth = FALSE){
    qs0 <- quantilesLogConDens(1 - x, res0)
    
    if (identical(smooth, FALSE)){c <- 3}
    if (identical(smooth, TRUE)){c <- 5}

    F1 <- evaluateLogConDens(qs0[, "quantile"], res1, which = c, gam = NULL, print = FALSE)

    if (identical(smooth, FALSE)){roc.x <- 1 - F1[, "CDF"]}
    if (identical(smooth, TRUE)){roc.x <- 1 - F1[, "smooth.CDF"]}

    return(as.numeric(roc.x))
}
