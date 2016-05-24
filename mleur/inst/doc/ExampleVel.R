library(lattice)
xyplot(vel, lwd=1.5, type="b", cex=0.7, pch=16, aspect=0.8,
 xlab="year", ylab="money velocity")
mleurDiag(vel)
mleur(vel)
dftest(vel)
