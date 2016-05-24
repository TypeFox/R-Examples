llmult <-
function (p, g) 
{
    index <- cbind(seq(along = g), as.numeric(g))
    p <- p[index]
    -2 * sum(log(p[p > .Machine$double.eps]))
}

