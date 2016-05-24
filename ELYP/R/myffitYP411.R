myffitYP411 <- function(x1, myY, myd, myZ, beta2) {
         tempT <- fitYP41(Y=myY, d=myd, Z=myZ, beta1=x1, beta2=beta2)
         return( - tempT$EmpLik)
}


