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
# ptmtotal <- proc.time()
cat("\n")
F2 <- planor.factors( factors=c(LETTERS[1:4], "bloc"),
                     nlevels=c(6,6,4,2,6),
                     block=~bloc)
M2 <- planor.model( model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D )
# ptm <- proc.time()
K2 <- planor.designkey(factors=F2, model=M2, nunits=144,
                       base=~A+B+D, max.sol=2)
cat("\n")
# ptm <- proc.time()
P2 <- planor.design(key=K2, select=c(1,1))
# ptm <- proc.time()
resum <- summary(K2[1,1])
#cat("Etude des alias\n")
#ptm <- proc.time()
#alias.designkey(K2[c(1,1)], model=M2[[1]][[1]])
#cat("TEMPS alias.designkey(K2[c(1,1)], model=M2[[1]][[1]])", proc.time()-ptm,"\n")
# cat("TEMPS summary(pick)", proc.time()-ptm,"\n")
# cat("TEMPS total", proc.time()-ptmtotal,"\n")
# Idem que
resum.l <- summary(K2)
                                        # REMARK: The following lines also work; they illustrate that the basic factors
# need not be part of the model but they must have been declared in planor.factors:
#
# F1 <- planor.factors( factors=c(LETTERS[1:4], "bloc", "BASE"), nlevels=rep(3,6) )
# M1 <- planor.model( model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D )
# K1 <- planor.designkey(factors=F1, model=M1, nunits=3^3,
#                       base=c("A","B","BASE"), max.sol=2)
#
