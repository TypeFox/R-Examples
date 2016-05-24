library(flsa)
library(testthat)

## tests using connlist objects
## just a repeat of a two-dimensional situation

connList <- vector("list", 4)
y <- 1:4

class(connList) = "connListObj"
connList[[1]] = as.integer(c(1,2))
connList[[2]] = as.integer(c(0,3))
connList[[3]] = as.integer(c(3,0))
connList[[4]] = as.integer(c(2,1))
names(connList) <- as.character(0:3) ## not necessary, just for illustration

res <- flsa(y, connListObj=connList)
res2 <- flsa(matrix(y, nrow=2))

expect_equal(res$BeginLambda, res2$BeginLambda, threshold=0.00001)
expect_equal(as.vector(flsaGetSolution(res, lambda2=c(0, 0.5, 1))),
             as.vector(flsaGetSolution(res2, lambda2=c(0, 0.5, 1))))

