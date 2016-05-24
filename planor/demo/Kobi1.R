library("planor")
#---------------------------------------------------------------------------
# EXAMPLES FROM THE PLANOR MANUAL
#---------------------------------------------------------------------------
# Exemple 1 page 7
#---------------------------------------------------------------------------
cat("\n")
cat("***************** EXEMPLE 1 PAGE 7 *****************\n")
cat("\n")
cat("Four 3-level treatment factors and one 3-level block factor\n")
cat("Model: bloc+(A+B+C+D)^2\n")
cat("Estimate: A+B+C+D\n")
cat("N=3^3\n")
cat("\n")
# 
cat("*** RUN ***\n")
F1 <- planor.factors( factors=c(LETTERS[1:4], "bloc"), nlevels=rep(3,5) )
M1 <- planor.model( model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D )
K1 <- planor.designkey(factors=F1, model=M1, nunits=3^3,
                       base=~A+B+C, max.sol=2, verbose=TRUE)
P1 <- planor.design(key=K1, select=2)
resum <- summary(K1[1])
