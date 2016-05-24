med2way <- function(formula, data){
#
#  Perform a J by K (two-way) anova on  medians where
#  all jk groups are independent.
#

  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  J <- nlevels(mf[,2])
  K <- nlevels(mf[,3])
  p <- J*K
  grp <- c(1:p)
  alpha <- .05
  
  x <- tapply(mf[,1], list(mf[,2], mf[,3]), function(xx) list(xx))
  if(p!=length(x)){
    print("Warning: The number of groups in your data is not equal to JK")
  }
  xbar<-0
  h<-0
  d<-0
  R<-0
  W<-0
  d<-0
  r<-0
  w<-0
  nuhat<-0
  omegahat<-0
  DROW<-0
  DCOL<-0
  xtil<-matrix(0,J,K)
  aval<-matrix(0,J,K)
  for (j in 1:p){
    #if(sum(duplicated(x[[grp[j]]]))>0)print("WARNING: TIED VALUES")
    xbar[j]<-median(x[[grp[j]]])
    h[j]<-length(x[[grp[j]]])
    d[j]<-msmedse(x[[grp[j]]], sewarn=FALSE)^2
  }
  d<-matrix(d,J,K,byrow=T)
  xbar<-matrix(xbar,J,K,byrow=T)
  h<-matrix(h,J,K,byrow=T)
  for(j in 1:J){
    R[j]<-sum(xbar[j,])
    nuhat[j]<-(sum(d[j,]))^2/sum(d[j,]^2/(h[j,]-1))
    r[j]<-1/sum(d[j,])
    DROW[j]<-sum(1/d[j,])
  }
  for(k in 1:K){
    W[k]<-sum(xbar[,k])
    omegahat[k]<-(sum(d[,k]))^2/sum(d[,k]^2/(h[,k]-1))
    w[k]<-1/sum(d[,k])
    DCOL[k]<-sum(1/d[,k])
  }
  D<-1/d
  for(j in 1:J){
    for(k in 1:K){
      xtil[j,k]<-sum(D[,k]*xbar[,k]/DCOL[k])+sum(D[j,]*xbar[j,]/DROW[j])-
        sum(D*xbar/sum(D))
      aval[j,k]<-(1-D[j,k]*(1/sum(D[j,])+1/sum(D[,k])-1/sum(D)))^2/(h[j,k]-3)
    }
  }
  Rhat<-sum(r*R)/sum(r)
  What<-sum(w*W)/sum(w)
  Ba<-sum((1-r/sum(r))^2/nuhat)
  Bb<-sum((1-w/sum(w))^2/omegahat)
  Va<-sum(r*(R-Rhat)^2)/((J-1)*(1+2*(J-2)*Ba/(J^2-1)))
  Vb<-sum(w*(W-What)^2)/((K-1)*(1+2*(K-2)*Bb/(K^2-1)))
  sig.A<-1-pf(Va,J-1,9999999)
  sig.B<-1-pf(Vb,K-1,9999999)
  # Next, do test for interactions
  Vab<-sum(D*(xbar-xtil)^2)
  dfinter<-(J-1)*(K-1)
  sig.AB<-1-pchisq(Vab,dfinter)
  #list(test.A=Va,p.val.A=sig.A,test.B=Vb,p.val.B=sig.B,test.AB=Vab,p.val.AB=sig.AB)
  result <- list(Qa=Va, A.p.value=sig.A, Qb=Vb, B.p.value=sig.B, Qab=Vab, AB.p.value=sig.AB, call = cl, varnames = colnames(mf))
  class(result) <- c("t2way")
  result
}
