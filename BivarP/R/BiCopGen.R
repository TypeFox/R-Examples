BiCopGen <- function(x, rodiny=c("weibull","weibull"),rodina="gumbel",
                     No, cens=TRUE, bicens=TRUE, digi=0)
{
 # x is a vector length 5. x[1:2] are parameters of the first marginal distribution. x[3:4] are parameters of the second one. x[5] is a copula parameter.
  DdX <- rep(1,No)
  DdY <- rep(1,No)
  a = x[length(x)]
  kopule <- archmCopula(rodina, param=a, dim=2)
  if(rodiny[1]=="weibull" | rodiny[1]=="gamma") {
    Lst1 = list(shape=x[1],scale=x[2])
  } else {if(rodiny[1]=="norm") {Lst1 = list(mean=x[1],sd=x[2])
  } else {if(rodiny[1]=="lnorm") {Lst1 = list(meanlog=x[1],sdlog=x[2])}}}
  if(rodiny[2]=="weibull" | rodiny[2]=="gamma") {
    Lst2 = list(shape=x[3],scale=x[4])
  } else {if(rodiny[2]=="norm") {Lst2 = list(mean=x[3],sd=x[4])
  } else {if(rodiny[2]=="lnorm") {Lst2 = list(meanlog=x[3],sdlog=x[4])}}}
  Lst = list(Lst1,Lst2) # podle rodiny
  bikop <- mvdc(kopule,rodiny,Lst)
  if (digi==0) XY <- ceiling(rMvdc(No,bikop)) else XY <- round(rMvdc(No,bikop),digi)
  if(cens==TRUE){
    bicop <- mvdc(kopule,rodiny,Lst)
    if (digi==0) CXY <- ceiling(rMvdc(No,bicop)) else CXY <- round(rMvdc(No,bicop),digi)
    if(bicens==TRUE){
      ixy <- which((CXY[,1] < XY[,1]) & (CXY[,2] < XY[,2]))
      if(length(ixy)>0) {
        XY[,1][ixy] <- CXY[,1][ixy]
        XY[,2][ixy] <- CXY[,2][ixy]
        DdX[ixy] <- 0
        DdY[ixy] <- 0
      }
    } else {
      ixc <- which(CXY[,1] < XY[,1])
      iyc <- which(CXY[,2] < XY[,2])
      if(length(ixc)>0) { XY[,1][ixc] <- CXY[,1][ixc]; DdX[ixc] <- 0 }
      if(length(iyc)>0) { XY[,2][iyc] <- CXY[,2][iyc]; DdY[iyc] <- 0 }
    }
  }
 return(list(X=XY[,1],Y=XY[,2],dX=DdX, dY=DdY, mvdc=bikop))
}