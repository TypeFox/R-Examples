isfe.fts <- function (data, max.order = N - 3, N = 10, h = 5:10, method = c("classical", 
    "M", "rapca"), mean = TRUE, level = FALSE, fmethod = c("arima", "ar",
    "ets", "ets.na", "struct", "rwdrift", "rw", "arfima"), lambda = 3, 
    ...) 
{
    method <- match.arg(method)
    fmethod <- match.arg(fmethod)
    isfe.vec <- matrix(NA, max.order + 1, length(h))
    for (i in 0:max.order) {
        isfe.vec[i + 1, ] <- colMeans(aveISFE(data, order = i, 
            N = N, h = h, method = method, fmethod = fmethod, 
            lambda = lambda, mean = mean, level = level, ...))
    }
    rownames(isfe.vec) <- 0:max.order
    colnames(isfe.vec) <- h
    return(isfe.vec)
}
