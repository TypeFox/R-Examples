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
# ptmtotal <- proc.time()
cat("\n")

#Rprof("Kobi0.designkey.Rprof.out")
# ptm <- proc.time()
K0 <- planor.designkey(factors=c(LETTERS[1:4], "bloc"),
                       block=~bloc,
                       nlevels=rep(2,5),
                       model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D,
                       nunits=2^3,
                       base=~A+B+C, max.sol=2, verbose=T)
# cat("TEMPS planor.designkey", proc.time()-ptm,"\n")
cat("\n")
#Rprof(NULL);Rprof("Kobi0.design.Rprof.out")
# ptm <- proc.time()
P0 <- planor.design(key=K0, select=2)
# cat("TEMPS design", proc.time()-ptm,"\n")
#Rprof(NULL);Rprof("Kobi0.summary.Rprof.out")
# ptm <- proc.time()
resum <- summary.designkey(K0[1])
# cat("TEMPS summary", proc.time()-ptm,"\n")
# cat("TEMPS total", proc.time()-ptmtotal,"\n")


# REMARK: The following lines also work; they illustrate that the basic factors
# need not be part of the model but they must have been declared in planor.factors:
#
F0 <- planor.factors( factors=c(LETTERS[1:4], "bloc", "BASE"), nlevels=rep(3,6),
                     block=~bloc )
K0 <- planor.designkey(factors=F0,
                        model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D,
                       nunits=3^3,   base=~A+B+BASE, max.sol=2)

