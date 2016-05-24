library("planor")
#---------------------------------------------------------------------------
# EXAMPLES FROM THE PLANOR MANUAL
#---------------------------------------------------------------------------
# Exemple 2 page 12
#---------------------------------------------------------------------------
cat("\n")
cat("***************** EXEMPLE 2 PAGE 12 *****************\n")
cat("\n")
cat("Four treatment factors at 6, 6, 4, 2 levels and one 6-level block factor\n")
cat("Model: bloc+(A+B+C+D)^2\n")
cat("Estimate: A+B+C+D\n")
cat("N=144=2^4 x 3^2\n")
cat("\n")
# 
cat("*** RUN ***\n")
F2 <- planor.factors( factors=c(LETTERS[1:4], "bloc"), nlevels=c(6,6,4,2,6) )
M2 <- planor.model( model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D )
K2 <- planor.designkey(factors=F2, model=M2, nunits=144,
                       base=~A+B+D, max.sol=2)
P2 <- planor.design(key=K2, select=c(1,1))
resum <- summary(K2[1,1])
# Same as:
resum.l <- summary(K2, selection=c(1,1))
