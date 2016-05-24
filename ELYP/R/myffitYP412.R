myffitYP412 <- function(x2, myY, myd, myZ, beta1) {
         tempT <- fitYP41(Y=myY, d=myd, Z=myZ, beta1=beta1, beta2=x2)
         return( - tempT$EmpLik)
}
 