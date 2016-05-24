#Author: Arni Magnusson, normally distributed in gplots library
richColors <- function (n, alpha = 1, rgb = FALSE) 
{
    if (n <= 0) 
        return(character(0))
    x <- seq(0, 1, length = n)
	r <- 1/(1 + exp(20 - 35 * x))
	g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
	b <- dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))

    rgb.m <- matrix(c(r, g, b), ncol = 3, dimnames = list(NULL, 
        c("red", "green", "blue")))
    col <- mapply(rgb, r, g, b, alpha)
    if (rgb) 
        attr(col, "rgb") <- cbind(rgb.m, alpha)
    return(col)
}

