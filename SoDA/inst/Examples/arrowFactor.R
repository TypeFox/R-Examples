arrowColors<- function(v, colors) {
    n <- length(v); k <- length(colors)
    midPoints <- (v[-1] + v[-n])/2
    colors[cut(midPoints, k)]
}
