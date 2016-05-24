library("planor")

cat("\n")
cat("***************** PLAN DE L'ARTICLE JTB09 *****************\n")
cat("********** Lurette, Touzeau, Lamboni, Monod ***********\n")
cat("\n")
cat("Eighteen 4-level treatment factors\n")
cat("N=2^12\n")
# 
cat("*** RUN ***\n")
# ptmtotal <- proc.time()
# Rprof("ST2.prof")
cat("\n")


cat("resolution 2\n")
cat("\n")
#  ptm <- proc.time()
ST.K <- planor.designkey(factors=LETTERS[1:18], nlevels=2,
                         model=~(A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R)^2 ,
                         nunits=2^12,
                       base=~A+B+C+D+E+F+G+H+I+J+K+L, max.sol=2)
# cat("*** TEMPS planor.designkey", proc.time()-ptm,"\n")
cat("\n")
ST.P <- planor.design(key=ST.K, select=2)
# cat("*** TEMPS design", proc.time()-ptm,"\n")
st2<- pick(ST.K,1)
# Increase print amount:
options(planor.max.print=100)
# ptm <- proc.time()
summary(st2)
# cat("*** TEMPS summary(pick)", proc.time()-ptm,"\n")
# cat("*** TEMPS total", proc.time()-ptmtotal,"\n")
#Rprof(NULL)

