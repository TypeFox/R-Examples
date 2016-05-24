myffitYP41 <- function(x, myY, myd, myZ) {
     x1 <- x[1]
     x2 <- x[2]
     tempT <- fitYP41(Y=myY, d=myd, Z=myZ, beta1=x1, beta2=x2)
     return( - tempT$EmpLik)
} 