"nsibs" <-
function(ls,lo,alpha, beta)
{
# calculate number of affected sibpairs needed to detect linkage to a 
# susceptibility gene attributed with a user specified lambda_sib/
# lambda_1.  Lambda_sib is the recurrence risk ratio for full- sibs
# (includes dominance variance), lambda_1 is the recurrence risk ratio
# for parent- offspring pairs (so does not include dominance variance).
# In these calculations we assume that parents and sibs are completely
# genotyped with markers that perfectly define the IBD configurations.
# alpha and beta correspond to the type 1 and 1 - type 2 error rates
# respectively.
  z0 <- 0.25/ls
  z1 <- 0.50*lo/ls
  zalpha <- qnorm(1-alpha)
  zbeta <- qnorm(beta)
  y <- 1 - 0.5*z1 - z0
  sigma <- sqrt(4*y*(1-y))
  mu <- 2*y-1
  n <- ((zalpha - sigma*zbeta)^2)/(2*mu^2)
  return(n)
}

