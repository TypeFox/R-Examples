adaptmh<-function(tab1,tab2,Gamma=1,alpha=0.05,double=FALSE,inc=0.25){
  stopifnot((length(Gamma)==1)&(Gamma>=1))
  stopifnot((length(alpha)==1)&(alpha>0)&(alpha<1))
  stopifnot((length(inc)==1)&(inc>0))
  stopifnot((double==TRUE)|(double==FALSE))
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

  adaptmhInternal<-function(tab1,tab2,gamma=gamma,alpha=alpha,double=double){
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
    if (double){
      lg1<-length(g1)
      gd<-c(g1[1],as.vector(rbind(rep(0,lg1-1),g1[2:lg1])))
      gboth<-gconv(gd,g2)
    }
    else {gboth<-gconv(g1,g2)}
    val1<-(0:(length(g1)-1))
    val2<-(0:(length(g2)-1))
    valboth<-(0:(length(gboth)-1))
    names(gboth)<-valboth
    names(g1)<-val1
    names(g2)<-val2
    actuala<-sum(tab1[1,1,])
    if (double) actualb<-(2*actuala)+sum(tab2[1,1,])
    else actualb<-actuala+sum(tab2[1,1,])

    jt<-outer(g1,g2,"*")
    va<-outer(val1,rep(1,length(g2)),"*")
    if (double){vb<-outer(2*val1,val2,"+")}
    else{vb<-outer(val1,val2,"+")}

    c1<-rev(cumsum(rev(g1)))
    cboth<-rev(cumsum(rev(gboth)))

    onealpha<-function(alpha1){
      if (sum(c1<=alpha1)>=1) mina<-min(val1[c1<=alpha1]) else mina<-Inf
      if (sum(cboth<=alpha1)>=1) minb<-min(valboth[cboth<=alpha1]) else minb<-Inf
      if (sum(c1<=(alpha1/3))>=1) maxa<-min(val1[c1<=(alpha1/3)]) else maxa<-max(val1)
      if (sum(cboth<=(alpha1/3))>=1) maxb<-min(valboth[cboth<=(alpha1/3)]) else maxb<-max(valboth)

      obj<-function(a,b){
        peither<-sum(jt[(va>=a)|(vb>=b)])
        pb<-sum(gboth[valboth>=b])
        pa<-sum(g1[val1>=a])
        adif<-abs(pa-pb)
        o<-c(a,b,peither,pa,pb,adif)
        names(o)<-c("a","b","peither","pa","pb","adif")
        o
      }
      o1<-NULL
      if ((mina==Inf)&(minb==Inf)){
        warning("At least one of your tables has such small counts that an alpha-level test at this Gamma is not informative.")
        stop("Results are not significant at level alpha, but no results could be with these data and Gamma.")
      }
      else if ((mina<Inf)&(minb==Inf)) {
        maxb<-Inf
        for (a in mina:maxa){
          o1<-rbind(o1,obj(a,minb))
        }
      }
      else if ((mina==Inf)&(minb<Inf)) {
        maxa<-Inf
        for (b in minb:maxb){
          o1<-rbind(o1,obj(mina,b))
        }
      }
      else {
        for (a in mina:maxa){
          for (b in minb:maxb){
            o1<-rbind(o1,obj(a,b))
          }
        }
      }

      o1<-as.data.frame(o1)
      o2<-o1[o1$peither<=alpha1,]
      o2<-o2[order(o2$a,o2$b),]
      ua<-sort(unique(o2$a))
      o3<-NULL
      bestb<-Inf
      for (i in 1:length(ua)){
        oi<-o2[o2$a==ua[i],]
        wh<-which.min(oi$b)
        currentb<-oi$b[wh]
        if (currentb<bestb) {
          o3<-rbind(o3,oi[wh,])
          bestb<-currentb
        }
      }
      o4<-o3[which.min(o3$adif),]
      o4
    }
    o<-NULL
    for (i in 1:length(alpha)){
      o<-rbind(o,onealpha(alpha[i]))
    }
    result<-rep("Accept at",length(alpha))
    result[(actuala>=o$a)|(actualb>=o$b)]<-"Reject at"
    o<-cbind(o,result,alpha)
    rownames(o)<-alpha
    list(A=actuala,B=actualb,result=o)

  }


  ninc<-floor((Gamma-1)/inc)
  if (ninc==0) gammas<-1
  else gammas<-c(1,1+(inc*(1:ninc)))
  if (!is.element(Gamma,gammas)) gammas<-c(gammas,Gamma)
  g1<-adaptmhInternal(tab1,tab2,gamma=gammas[1],alpha=alpha,double=double)
  out<-g1$result
  acta<-g1$A
  actb<-g1$B
  if (length(gammas)>=2) {
    for (i in 2:length(gammas)) {
      out<-rbind(out,adaptmhInternal(tab1,tab2,gamma=gammas[i],alpha=alpha,double=double)$result)
    }
  }
  if (dim(out)[1]==1) out<-as.matrix(out,1,length(out))
  out<-as.data.frame(out)
  rownames(out)<-gammas
  Gamma<-gammas
  out<-cbind(out,Gamma)
  if (!all(out$result=="Reject at")) {
    first<-which(out$result=="Accept at")[1]
    out<-out[1:first,]
  }
  list(A=acta,B=actb,result=out)
}
