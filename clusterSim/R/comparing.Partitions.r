comparing.Partitions<-function(cl1,cl2,type="nowak"){
  if(!is.vector(cl1) || !is.vector(cl2)){
    stop("cl1 or cl2 should be vectors")
  }
  if(length(cl1)!=length(cl2)){
    stop("cl1 and cl2 should have equal sizes")
  }
  if(!any(type==c("rand","crand","nowak"))){
    stop("type paramater should be one the following: rand, crand, nowak")
  }
  if(type=="rand"){
    resul<-classAgreement(table(cl1,cl2))$rand
  }
  if(type=="crand"){
    resul<-classAgreement(table(cl1,cl2))$crand
  }
  if(type=="nowak"){
    ktab<-tab<-table(cl1,cl2)
    for(i in 1:nrow(tab)){
    for(j in 1:ncol(tab)){
      ktab[i,j]<-tab[i,j]/max(sum(tab[i,]),sum(tab[,j]))
    }
    }
    resul<-sum(apply(ktab,1,max),apply(ktab,2,max))/(length(unique(cl1))+length(unique(cl2)))
  }
  resul
}
#print(.comparing.Partitions(c("b","c","b","c","b","a"),c("b","a","a","b","c","c")),type="nowak")
