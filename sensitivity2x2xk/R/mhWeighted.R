mhWeighted<-function(tab1,tab2,Gamma=1){
  stopifnot((length(Gamma)==1)&(Gamma>=1))
  d1<-dim(tab1)
  d2<-dim(tab2)
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
      mx<-min(n,m1)
      mn<-max(0,n-m2)
      x<-mn:mx
      g<-rep(0,mx+1)
      pr<-dFNCHypergeo(x,m1,m2,n,gamma)
      g[(mn+1):(mx+1)]<-pr
      g
    }

    one2x2xk<-function(tbk,gamma=gamma){
      k<-dim(tbk)[3]
      g<-1
      for (i in 1:k){
        gi<-one2x2(tbk[,,i],gamma=gamma)
        g<-gconv(g,gi)
      }
      g
    }

    g1=one2x2xk(tab1,gamma=gamma)
    g2=one2x2xk(tab2,gamma=gamma)
    lg1<-length(g1)
    gd<-c(g1[1],as.vector(rbind(rep(0,lg1-1),g1[2:lg1])))
    gboth<-gconv(gd,g2)
    val1<-(0:(length(g1)-1))
    val2<-(0:(length(g2)-1))
    valboth<-(0:(length(gboth)-1))
    names(gboth)<-valboth
    names(g1)<-val1
    names(g2)<-val2
    actuala<-sum(tab1[1,1,])
    actualb<-(2*actuala)+sum(tab2[1,1,])
    pval<-sum(gboth[valboth>=actualb])
    list(pval=pval,WeightedMH=actualb)
}
