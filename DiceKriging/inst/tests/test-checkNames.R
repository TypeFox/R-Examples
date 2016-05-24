c1 <- rep(1, 5)
c2 <- rep(2, 5)
c3 <- rep(3, 5)

X1 <- data.frame(x1=c1, x2=c2, x3=c3)

permutation.list <- list(c(1, 2, 3), c(1, 3, 2), c(3, 2, 1), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2))
precision <- 1e-8

for (i in 1:6) {
  permutation <- permutation.list[[i]]
  permutation.name <- paste(permutation, collapse="")
  X2 <- X1[, permutation]  # Does also permute the names
  newdata <- checkNames(X1=X1, X2=X2)
  departure_newdata_X1 <- sum((X1 - newdata)^2)
  test_that(desc=paste("test checkNames, permutation=", permutation.name), 
            expect_true(departure_newdata_X1 < precision)) 
}

