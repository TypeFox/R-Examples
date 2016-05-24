magpieComp<-function(bench,comp,reg=NA) {

  if(!is.na(reg)){
    bench<-bench[reg,,]
    comp<-comp[reg,,]
  }
  # keep only common variable names
  comp1<-collapseNames(comp[,,intersect(getNames(comp),getNames(bench))])
  bench1<-collapseNames(bench[,,intersect(getNames(comp),getNames(bench))])
  
  nas<-which(is.na(comp1)&is.na(bench1))
  bench1[nas]<-0
  comp1[nas]<-0
  rem0<-which(bench1==0 & comp1==0)
  bench1[rem0]<-1
  comp1[rem0]<-1
  data<-200*abs(bench1-comp1)/(abs(bench1)+abs(comp1))
  
  out<-list()
  out$varnames<-list()
  out$varnames$bench_only<-setdiff(getNames(bench),getNames(comp))
  out$varnames$comp_only<-setdiff(getNames(comp),getNames(bench))
  d<-as.data.frame(data,rev=2)
  names(d)<-sub(".value","value",names(d))
  out$diff_sorted<-d[order(d$value,decreasing = TRUE,na.last = T),]
  out$data<-data
  
return(out)
}