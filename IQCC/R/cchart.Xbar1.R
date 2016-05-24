cchart.Xbar1 <- function(x, sizes)
{
    a <- rowMeans(x)
    x2bar <- mean(a)
    sigma <- sd.xbar(x)
    qcc(x, type = "xbar", sizes)
	stat <- list(c(x2bar, sigma))
    return(stat)
}