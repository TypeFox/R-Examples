#library(r2stl)

# Let's do the classic persp() demo plot

x <- seq(-10, 10, length= 100)

y <- x

f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }

z <- outer(x, y, f)

z[is.na(z)] <- 1

r2stl(x, y, z, filename="lovelyfunction.stl", show.persp=TRUE)

