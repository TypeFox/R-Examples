##
##  c r o s s n . R  tests
##

crossn <- pracma::crossn

x <- c(1.0, 0.0, 0.0)
y <- c(1.0, 0.5, 0.0)
z <- c(0.0, 0.0, 1.0)
identical(pracma::dot(x, crossn(rbind(y, z))),
          det(rbind(x, y, z)))
