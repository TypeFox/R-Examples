# Tukey HSD method applied to Kruskal-Wallis test.
# resp is the response variable, and grp indicates group.
tukey.kruskal.test<-function(resp,grp,alpha=.05){
  ld<-split(resp,grp)
  nv<-sapply(ld,length)
  rm<-sapply(split(rank(resp),grp),mean)
  ng<-length(rm)
  N<-sum(nv)
  out<-array(NA,c(ng*(ng-1)/2,2))
  dimnames(out)<-list(rep("",ng*(ng-1)/2),NULL)
  rr<-rank(resp)
  count<-0
  # This is similar to higgins.fisher.kruskal.test, except that
  # here a loop is substituted for vector calc.,
  # and we calculate confidence intervals instead of p-values.
  for(i in 1:(ng-1)) for(j in (i+1):ng){
    count<-count+1
    out[count,]<-rm[j]-rm[i]+c(-1,1)*qtukey(1-alpha,ng,N-ng)*
      sqrt(N*(N+1)/24)*sqrt(1/nv[i]+1/nv[j])
    dimnames(out)[[1]][count]<-paste(i,j,sep="-")
  }
  different<-(out[,1]>0)|(out[,2]<0)
  names(different)<-dimnames(out)[[1]]
  return(different)
}
