LRNormal2Mean <-
function(x,y,alpha){
  xbar <- mean(x); ybar <- mean(y)
  nx <- length(x); ny <- length(y)
  Sx <- var(x); Sy <- var(y)
  Sp <- ((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2)
  tcalc <- abs(xbar-ybar)/sqrt(Sp*(1/nx+1/ny))
  conclusion <- ifelse(tcalc>qt(df=nx+ny-2,p=alpha/2),
                       "Reject Hypothesis H","Fail to Reject Hypothesis H")
  return(c(tcalc,conclusion,Sp))
}
