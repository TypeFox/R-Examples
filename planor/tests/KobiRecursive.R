library("planor")
#---------------------------------------------------------------------------
# EXAMPLE FROM "Automatic generation of asymmetrical regular designs"
#---------------------------------------------------------------------------
# Exemple 3 page 13
#---------------------------------------------------------------------------
cat("\n")
cat("***************** EXEMPLE 3 PAGE 13 of algo report *****************\n")
cat("\n")
cat("Row-Column (3x2) design with 2 units per row-column combination\n")
cat("3 treatment factors A(3) row-constant, B1 (2), B2 (2)\n")
cat("Model: R*C + (A+B1+B2)^2\n")
cat("Estimate: A:B1+A:B2\n")
cat("N=3x2x2\n")
cat("\n")
#
cat("*** RUN ***\n")
# ptmtotal <- proc.time()
cat("\n")
F2 <- planor.factors( factors=c("R","C","U","A","B1","B2"), nlevels=c(3,2,2,3,2,2) )
M2 <- planor.model( model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2 )
# ptm <- proc.time()
K2 <- planor.designkey(factors=F2, model=M2, nunits=12,
                       base=~R+C+U, max.sol=2)
# cat("TEMPS planor.designkey", proc.time()-ptm,"\n")
cat("\n")
# ptm <- proc.time()
P2 <- planor.design(K2)
# cat("TEMPS design", proc.time()-ptm,"\n")
resum <- summary(K2)

# cat("TEMPS total", proc.time()-ptmtotal,"\n")

