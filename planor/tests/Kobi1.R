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
# ptmtotal <- proc.time()
cat("\n")
F1 <- planor.factors( factors=c(LETTERS[1:4], "bloc"), nlevels=rep(3,5),
                     block=~bloc)
M1 <- planor.model( model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D )
# ptm <- proc.time()
K1 <- planor.designkey(factors=F1, model=M1, nunits=3^3,
                       base=~A+B+C, max.sol=2, verbose=TRUE)
# cat("TEMPS planor.designkey", proc.time()-ptm,"\n")
cat("\n")
# ptm <- proc.time()
P1 <- planor.design(key=K1, select=2)
# cat("TEMPS design", proc.time()-ptm,"\n")
# ptm <- proc.time()
resum <- summary(K1[1])
# cat("TEMPS summary(pick)", proc.time()-ptm,"\n")
#cat("Etude des alias\n")
#ptm <- proc.time()
#alias.designkey(K1[1], model=M1[[1]][[1]])
#cat("TEMPS alias.designkey(K1[1], model=M1[[1]][[1]])", proc.time()-ptm,"\n")
# cat("TEMPS total", proc.time()-ptmtotal,"\n")


                                        # REMARK: The following lines also work; they illustrate that the basic factors
# need not be part of the model but they must have been declared in planor.factors:
#
# F1 <- planor.factors( factors=c(LETTERS[1:4], "bloc", "BASE"), nlevels=rep(3,6) )
# M1 <- planor.model( model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D )
# K1 <- planor.designkey(factors=F1, model=M1, nunits=3^3,
#                       base=c("A","B","BASE"), max.sol=2)
#
