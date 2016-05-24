# Fisher's LSD method applied to the Kruskal-Wallis test, as presented by Higgins.
higgins.fisher.kruskal.test<-function(resp,grp,alpha=.05){
  ld<-split(resp,grp)
  grpl<-unique(grp)
  nv<-sapply(ld,length)
  rm<-sapply(split(rank(resp),grp),mean)
  ng<-length(rm)
  N<-sum(nv)
  out<-NULL
  if(kruskal.test(resp,grp)$p.value<alpha){
    pvfs<-(1-pnorm(abs(outer(rm,rm,"-"))/
                     sqrt(N*(N+1)*outer(1/nv,1/nv,"+")/12)))*2+
      (outer(seq(ng),seq(ng),"-") <=0)
    i1<-outer(seq(ng),rep(1,ng))
    out<-cbind(t(i1)[pvfs<alpha],i1[pvfs<alpha])
  }
  return(out)
}
