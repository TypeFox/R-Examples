test.create.RandomMatrixEstimator <- function()
{
  a <- create(RandomMatrixEst, hint=c(6,2))
  checkTrue(isa(RandomMatrixEst,a))
  checkTrue(a$hint == c(6,2))
}


