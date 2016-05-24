myffitYP3 <- function(x, myY, myd, myZ) {
x1 <- x[1]
x2 <- x[2]
x3 <- x[3]
tempT <- fitYP3(Y=myY,d=myd,Z=myZ, beta1=c(x3, x1), beta2= c(x3, x2), lam=0, fun=function(t){as.numeric(t <= 0.175)} )
return(- tempT$LogEmpLik)
}
