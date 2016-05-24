# Julia Bischof
# 12-08-2015

#library(ape)

dist.PCoA<-function(dist.tab=NULL, correction=c("lingoes","cailliez","none")){
  if(length(correction)==0){
    correction=="none"
  }
  out<-list()
  if(is.list(dist.tab)==T){
    for(i in 1:length(dist.tab)){
      out<-c(out,list(tryCatch(pcoa(as.dist(dist.tab[[i]]),correction=correction),
                               error=function(e){"no results available"})))
    }
    names(out)<-names(dist.tab)
  }else{
    out<-list(pcoa(as.dist(dist.tab),correction=correction))
  }
  return(out)
}

