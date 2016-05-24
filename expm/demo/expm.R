require("expm")


#diagonalisable matrix
T <- cbind(c(-2, -3, 2), c(2, -2, 1), c(0, 2, -2))

expm(T)


# numerically singular matrix
T <- cbind(c(-2, 0, 0), c(2, -2, 0), c(0, 2, -2))
expm(T)

#solve shows T is numerically singular
try(solve(eigen(T)$vectors))

#singular matrix
T <- cbind(c(0, 0, 0), c(2, 0, 0), c(1, 2, 0))
expm(T)
stopifnot(all.equal(logm(expm(T)), T))
## and show how close it is
all.equal(logm(expm(T)), T,  tolerance=0)# 2.39e-15 {64b ubuntu 12-04}
