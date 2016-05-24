## assumes convex objectives
`profSmooth` <-
function (prof, ...) 
UseMethod("profSmooth")
`profSmooth.profileModel` <-
function (prof, n.interpolations = 100, ...) 
{
    isNA <- prof$isNA
    profRes <- prof$profiles
    p <- length(profRes)
    BetasNames <- names(profRes)
    intersects <- prof$intersects
    quantile <- prof$quantile
    result <- matrix(rep(c(-Inf, Inf), each = p), p, 2)
    for (i in 1:p) {
        if (isNA[i]) {
            result[i, ] <- NA
            next
        }
        profRes.i <- profRes[[i]]
        smoothed <- spline(profRes.i, n = n.interpolations)
        min.which <- which.min(smoothed$y)
        bb <- smoothed$x[min.which]
        left <- which(smoothed$x < bb)
        right <- which(smoothed$x >= bb)
        if (intersects[i, 1]) 
            result[i, 1] <- approx(x = smoothed$y[left], y = smoothed$x[left], 
                xout = quantile)$y
        if (intersects[i, 2]) 
            result[i, 2] <- approx(x = smoothed$y[right], y = smoothed$x[right], 
                xout = quantile)$y
    }
    dimnames(result) <- list(BetasNames, c("Lower", "Upper"))
    result
}
