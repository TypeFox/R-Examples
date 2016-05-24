input <- matrix(runif(6), 3, 2)
target <- rowSums(input)

darch <- new_dnn(c(2, 1))
predict(darch, newdata = input)

x <- cbind(input, 1)
weight <- getLayer(darch, 1)[[1]]

y <- x %*% weight

y
