kern.B <- function (x, xi, h, g = 0)
{
    ((1-((xi-x)/h)^2)^g)/beta(.5,g+1)/h
}
