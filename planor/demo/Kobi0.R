library("planor")
#---------------------------------------------------------------------------
# EXAMPLES FROM THE PLANOR MANUAL
#---------------------------------------------------------------------------
# Example 1 page 2
#---------------------------------------------------------------------------
cat("\n")
cat("***************** EXAMPLE 1 PAGE 2 *****************\n")
cat("\n")
cat("Four 2-level treatment factors and one 2-level block factor\n")
cat("Model: bloc+A+B+C+D\n")
cat("Estimate: A+B+C+D\n")
cat("N=2^3\n")
cat("\n")
# 
cat("*** RUN ***\n")
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"),
                       nlevels=rep(2,5),
                       model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
                       nunits=2^3,
                       base=~A+B+C, max.sol=2, verbose=T)
P0 <- planor.design(key=K0, select=2)
resum <- summary(K0[1])
