library(flsa)
library(testthat)

y <- as.numeric(1:10)

## run a simple example with flsa
res.usingLambdas <- flsa(y, lambda2=c(0, 1, 3, 6, 10, 12.5))
res.solObj <- flsa(y)
res.solFromObj <- flsaGetSolution(res.solObj, lambda2=c(0, 1, 3, 6, 10, 12.5))

## then with flsaTopDown
res.topdown <- flsaTopDown(y)

## and the same using integer input
res.integer <- flsa(1:10, lambda2=c(0, 1, 3, 6, 10, 12.5))


## check that the solutions are all the same and correct
## first check that the solution objects contain the correct lambda values
expect_equal(sort(unique(round(res.solObj$mergeLambda,12)), decreasing=FALSE), c(-1, 1, 3, 6, 10, 12.5))
expect_equal(sort(unique(round(res.topdown$Lambdas, 12)), decreasing=FALSE), c(0, 1, 3, 6, 10, 12.5))

topdown.matrix <- res.topdown$Solution[as.character(c(0,1,3,6,10,12.5)),]
colnames(topdown.matrix) <- NULL
rownames(topdown.matrix) <- NULL

expect_equal(max(topdown.matrix- res.usingLambdas), 0, tolerance=0.000001)
expect_equal(max(topdown.matrix- res.solFromObj), 0, tolerance=0.000001)
expect_equal(max(topdown.matrix- res.integer), 0, tolerance=0.000001)
