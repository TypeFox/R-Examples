##
##  i n t e r p 1 . R Test suite
##


interp1 <- pracma::interp1

x <- c(0.0, 0.5, 1.0, 1.5)
y <- x^2
xi <- c(0.25, 0.75, 1.25)

identical(interp1(x, y, xi, method="constant"), c(0.0, 0.25, 1.0))
identical(interp1(x, y, xi, method="linear"),   c(0.125, 0.625, 1.625))
identical(interp1(x, y, xi, method="nearest"),  c(0.25, 1.00, 2.25))
identical(interp1(x, y, xi, method="spline"),   c(0.0625, 0.5625, 1.5625))

# Not yet implemented
# identical(interp1(x, y, xi, method="cubic"),    c(0.0781, 0.5547, 1.5547))
