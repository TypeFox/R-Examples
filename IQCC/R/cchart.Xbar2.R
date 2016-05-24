cchart.Xbar2 <- function(x, x2bar, sigma, sizes)
{
    qcc(x, type = "xbar", center = x2bar, std.dev = sigma)
}