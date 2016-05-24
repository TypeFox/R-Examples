##
##  p i e c e w i s e . R  Test suite
##


piecewise <- pracma::piecewise

x <- c(0, 1,  2, 3, 4, 5)
y <- c(1, 1, -1, 0, 1, 0)
identical(piecewise(x, y)$area,  1.5)
identical(piecewise(x, y)$zeros, c(1.5, 3, 5))
identical(piecewise(x, y, abs = TRUE)$area,  3.0)
identical(piecewise(x, y, abs = TRUE)$zeros, c(1.5, 3, 5))
