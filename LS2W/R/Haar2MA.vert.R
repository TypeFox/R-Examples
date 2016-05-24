`Haar2MA.vert` <-
function(n, sd = 1., order = 5.)
{
#
# Generate realization of a vertical 2-D Haar MA field
#
# n = dimension of the square image
# sd = variance of orthonromal sequence
# order = MA order
#
temp <- rnorm(n = (n + (2.^order) - 1.)^2., mean = 0., sd = sd)
z <- matrix(temp, nrow = (n + (2.^order) - 1.))
x <- matrix(0., nrow = n, ncol = n)
J <- order
for(i in (2.^J):(2.^(J - 1.) + 1.)) {
   for(j in (2.^J):(2.^(J - 1.) + 1.))
      x <- x + z[i:(n + i - 1.), j:(n + j - 1.)]
}
for(i in (2.^J):(2.^(J - 1.) + 1.)) {
   for(j in (2.^(J - 1.)):1.)
      x <- x + z[i:(n + i - 1.), j:(n + j - 1.)]
}
for(i in (2.^(J - 1.)):1.) {
   for(j in (2.^J):(2.^(J - 1.) + 1.))
      x <- x - z[i:(n + i - 1.), j:(n + j - 1.)]
}
#       print(z)
for(i in (2.^(J - 1.)):1.) {
   for(j in (2.^(J - 1.)):1.)
      x <- x - z[i:(n + i - 1.), j:(n + j - 1.)]
}
x <- 2.^( - J) * x
return(x)
}

