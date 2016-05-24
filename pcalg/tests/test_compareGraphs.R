library(pcalg)

amat1 <- t(cbind(c(0,1,0,1,0),c(0,0,1,0,1),c(0,0,0,1,1),c(0,0,0,0,1),c(0,0,0,0,0)))
amat2 <- t(cbind(c(0,1,0,1,1),c(0,0,0,0,1),c(0,0,0,0,1),c(0,0,0,0,1),c(0,0,0,0,0)))

g1 <- as(amat1,"graphNEL")
g2 <- as(amat2,"graphNEL")

res <- compareGraphs(g1,g2)

if ((round(res["tpr"],5)!=0.83333) | (round(res["fpr"],5)!=0.5) | (round(res["tdr"],5)!=0.71429)) {
  stop("Test of compareGraphs: Theoretical values not matched!")
}

