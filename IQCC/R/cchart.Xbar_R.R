cchart.Xbar_R <- function(x, sizes)
{
    par(mfrow = c(1, 2))                 # setup 1 row and 2 columns for plotting
    qcc(x, type = "xbar", add.stats = F)
    qcc(x, type = "R", add.stat = F)
}