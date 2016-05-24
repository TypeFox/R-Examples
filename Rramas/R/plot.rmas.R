plot.rmas <-
function (x, sum = TRUE, mean = FALSE, type = "l", harvest=FALSE, ...) 
{
        if(length(names(x))>0){ # bifurcación provisional para separar rmas de projectn y de projectn2
	    if(harvest==FALSE) x<- x$vn else x <-x$harvest # seleciona trayectorai de la población o el manejo
        }
	 
    nl <- length(x)
    stages <- dim(x[[1]])[1]
    time <- dim(x[[1]])[2]
    if (nl == 1) {
        if (sum == TRUE) {
            plot(0:(time - 1), apply(x[[1]], 2, sum), type = type, 
                xlab = "time", ylab = "abundance", ...)
        }
        if (sum != TRUE) {
            matplot(0:(time - 1), t(x[[1]]), type = type, xlab = "time", 
                ylab = "abundance", col = rainbow(stages), pch = 1, 
                ...)
        }
    }
    if (nl > 1) {
        if (mean == TRUE) {
            if (sum == TRUE) {
                plot(0:(time - 1), apply(sapply(x, function(re) apply(re, 
                  2, sum)), 1, mean), type = type, xlab = "time", 
                  ylab = "abundance", ...)
            }
            if (sum != TRUE) {
                cosa <- NULL
                for (i in 1:stages) cosa <- rbind(cosa, apply(sapply(x, 
                  function(st) st[i, ]), 1, mean))
                matplot(0:(time - 1), t(cosa), type = type, xlab = "time", 
                  ylab = "abundance", col = rainbow(stages), 
                  pch = 1, ...)
            }
        }
        if (mean != TRUE) {
            if (sum == TRUE) {
                matplot(0:(time - 1), (sapply(x, function(re) apply(re, 
                  2, sum))), type = type, xlab = "time", ylab = "abundance", 
                  pch = 1, ...)
            }
            if (sum != TRUE) {
                par(mfrow = c(1, stages))
                for (i in 1:stages) matplot(0:(time - 1), sapply(x, 
                  function(st) st[i, ]), xlab = "time", ylab = "abundance", 
                  type = type, col = rainbow(stages)[i], pch = 1, 
                  main = paste("stage", i), ...)
            }
        }
    }
}
