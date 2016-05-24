sqndwdecomp <-
function (x, J, filter.number, family) 
{
    lx <- length(x)
    ans <- matrix(0, nrow = J, ncol = length(x))
    dw <- hwwn.dw(J, filter.number, family)
    longest.support <- length(dw[[J]])
    scale.shift <- 0
    for (j in 1:J) {
        l <- length(dw[[j]])
        init <- (filter.number - 1) * (lx - 2^j)
        for (k in 1:lx) {
            yix <- seq(from = k, by = 1, length = l)
            yix <- ((yix - 1)%%lx) + 1
            ans[j, k] <- sum(x[yix] * dw[[j]]^2)
        }
        if (filter.number == 1) 
            scale.shift <- 0
        else {
            scale.shift <- (filter.number - 1) * 2^j
        }
        ans[j, ] <- guyrot(ans[j, ], scale.shift)
    }
    return(ans)
}
