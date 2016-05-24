cutter <-
function (x, cuts) 
{
    ncuts <- length(cuts)
    xi <- cut(x, breaks = cuts, labels = FALSE, include.left = TRUE)
    xi[x <= cuts[1]] <- 0
    xi[x >= cuts[ncuts]] <- ncuts
    xi <- xi + 1
    invisible(xi)
}
