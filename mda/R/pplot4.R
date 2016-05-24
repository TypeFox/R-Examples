pplot4 <-
function (x, ...) 
{
    par(mfrow = c(3, 2))
    for (i in 1:3) {
        for (j in (i + 1):4) pplot(x[, c(i, j)], xlab = paste("var", 
            i), ylab = paste("var", j), ...)
    }
}

