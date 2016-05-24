#Difference in BAA and AAA corporate bonds 
library(lattice)
xyplot(DiffBA, lwd=1.5, type="b", cex=0.7, pch=16, aspect=0.8,
 xlab="year", ylab="money velocity")
mleurDiag(DiffBA)
mleur(DiffBA)
dftest(DiffBA)
ar1est(DiffBA)
ar1est(DiffBA, method="LSE")
