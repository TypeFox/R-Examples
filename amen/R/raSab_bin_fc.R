#' Simulate a and Sab from full conditional distributions under bin likelihood
#' 
#' Simulate a and Sab from full conditional distributions under bin likelihood
#' 
#' 
#' @usage raSab_bin_fc(Z, Y, a, b, Sab, SS = round(sqrt(nrow(Z))))
#' @param Z a square matrix, the current value of Z
#' @param Y square binary relational matrix
#' @param a current value of row effects
#' @param b current value of column effects
#' @param Sab current value of Cov(a,b)
#' @param SS number of iterations
#' @return \item{Z}{new value of Z} \item{Sab}{new value of Sab} \item{a}{new
#' value of a}
#' @author Peter Hoff
#' @export raSab_bin_fc
raSab_bin_fc <-
function(Z,Y,a,b,Sab,SS=round(sqrt(nrow(Z))))
{
  E<-Z-a%*%t(rep(1,nrow(Z))) 
  MEL<-MEU<- -E
  MEL[!is.na(Y) & Y==0]<- -Inf
  MEU[!is.na(Y) & Y==1]<-Inf 
  MEL[is.na(Y)]<- -Inf ; MEU[is.na(Y)]<- Inf
#  diag(MEU)<-Inf;diag(MEL)<- -Inf
  lba<-apply(MEL,1,max) 
  lba[is.na(lba)]<- -Inf
  uba<-apply(MEU,1,min) 
  uba[is.na(uba)]<- Inf

  for(ss in 1:SS)
  {
    ea<-b*Sab[1,2]/Sab[2,2]
    sa<-sqrt(Sab[1,1]-Sab[1,2]^2/Sab[2,2])
    a<-ea+sa*qnorm(runif(nrow(Z),pnorm((lba-ea)/sa),pnorm((uba-ea)/sa)))
    Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a,b))),3+nrow(Z)))
  }
  list(Z=E+a%*%t(rep(1,nrow(Z))),a=a,Sab=Sab)
}
