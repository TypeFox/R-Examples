gx.adjr2 <-
function (mr2, n, p) 
{
    adjr2 <- 1 - ((n - 1)/(n - p)) * (1 - mr2)
    adjr2 <- signif(adjr2, 3)
    cat(paste("  R^2 = ", mr2, ";  n = ", n, ";  p = ", p, sep = ""), 
        "\t--  Adjusted R^2 =", adjr2, "\n\n")
    invisible()
}
