myffitYP4 <- function(x, myY, myd, myZ) {
         if(length(x) != 2) stop("x must be a length 2 vector")
         if(dim(myZ)[2] == 1) stop("check dim of Zmat")
         x1 <- x[1]
         x2 <- x[2]
         tempT <- fitYP4(Y=myY, d=myd, Z=myZ, beta1=x1, beta2=x2)
         return( - tempT$EmpLik)
}
