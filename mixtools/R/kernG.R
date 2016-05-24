kern.G <- function (x, xi, h)
{
    exp(-((xi-x)/h)^2/2)/sqrt(2*pi)/h
}
