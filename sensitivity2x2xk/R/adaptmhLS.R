adaptmhLS<-function(tab1,tab2,Gamma=1,double=FALSE){
  #large sample version   Requires mvtnorm
  d1<-dim(tab1)
  d2<-dim(tab2)
  stopifnot(Gamma>=1)
  stopifnot((d1[1]==2)&(d1[2]==2)&(d2[1]==2)&(d2[2]==2))
  if (length(d1)==2) tab1<-array(tab1,c(2,2,1))
  if (length(d2)==2) tab2<-array(tab2,c(2,2,1))
  stopifnot(sum(as.vector(tab1))>0)
  stopifnot(sum(as.vector(tab2))>0)
  tab1<-tab1[,,apply(tab1,3,sum)>0]
  tab2<-tab2[,,apply(tab2,3,sum)>0]
  if (length(dim(tab1))==2) tab1<-array(tab1,c(2,2,1))
  if (length(dim(tab2))==2) tab2<-array(tab2,c(2,2,1))
  gamma<-Gamma

  check2x2xktable<-function(tab){
    r1<-tab[1,1,]+tab[1,2,]
    r2<-tab[2,1,]+tab[2,2,]
    c1<-tab[1,1,]+tab[2,1,]
    c2<-tab[1,2,]+tab[2,2,]
    mc<-pmin(pmin(r1,r2),pmin(c1,c2))
    if (max(mc)==0) {
      warning("One of the 2x2 or 2x2xk tables is degenerate.")
      stop("There is a problem with your data.")
    }
  }
  check2x2xktable(tab1)
  check2x2xktable(tab2)

  one2x2<-function(tb,gamma=gamma){
    m1<-tb[1,1]+tb[1,2]
    m2<-tb[2,1]+tb[2,2]
    n<-tb[1,1]+tb[2,1]
    expect<-meanFNCHypergeo(m1, m2, n, gamma)
    varian<-varFNCHypergeo(m1, m2, n, gamma)
    list(expect=expect,varian=varian)
  }

  one2x2xk<-function(tbk,gamma=gamma){
    k<-dim(tbk)[3]
    ex<-rep(NA,k)
    va<-rep(NA,k)
    for (i in 1:k){
      one<-one2x2(tbk[,,i],gamma=gamma)
      ex[i]<-one$expect
      va[i]<-one$varian
    }
    list(ex=ex,va=va)
  }

  mo1<-one2x2xk(tab1,gamma=gamma)
  ex1<-sum(mo1$ex)
  va1<-sum(mo1$va)
  mo2<-one2x2xk(tab2,gamma=gamma)
  ex2<-sum(mo2$ex)
  va2<-sum(mo2$va)
  if (double){
    exboth<-(ex1*2)+ex2
    vaboth<-(va1*4)+va2
    rho<-(2*va1)/sqrt(vaboth*va1)
    co<-matrix(c(1,rho,rho,1),2,2)
  }
  else {
    exboth<-ex1+ex2
    vaboth<-va1+va2
    rho<-va1/sqrt(vaboth*va1)
    co<-matrix(c(1,rho,rho,1),2,2)
    }
  actuala<-sum(tab1[1,1,])
  actual2<-sum(tab2[1,1,])
  if (double) actualb<-(2*actuala)+actual2
  else actualb<-actuala+actual2

  deva<-(actuala-ex1)/sqrt(va1)
  dev2<-(sum(actual2-ex2)/sqrt(va2))
  devboth<-(actualb-exboth)/sqrt(vaboth)
  maxdev<-max(deva,devboth)
  o1<-matrix(NA,3,5)
  rownames(o1)<-c("table1","table2","both")
  colnames(o1)<-c("statistic","expectation","variance","deviate","pvalue")
  o1[,1]<-c(actuala,actual2,actualb)
  o1[,2]<-c(ex1,ex2,exboth)
  o1[,3]<-c(va1,va2,vaboth)
  o1[,4]<-c(deva,dev2,devboth)
  o1[,5]<-1-pnorm(o1[,4])

  pval<-1-pmvnorm(lower=c(-Inf,-Inf),upper=c(maxdev,maxdev),corr=co)
  pval<-pval[1]
  list(pval=pval,maxdeviate=maxdev,correlation=co[1,2],detail=o1)

}
