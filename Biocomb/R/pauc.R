pauc <- function(auc,n=100,n.plus=0.5,labels=numeric(),pos=numeric()){

  pauc<-array(0,length(auc))
  if(length(pos)==0)
  {
  for(i in 1:length(auc)){
  np <- n.plus
  if (n.plus<1){
    if (n.plus<0){
      np <- round(0.5*n)
    }else{
      np <- round(n.plus*n)
    }
  }
  nm <- n - np
  pauc[i]<-pwilcox(nm*np*(1-auc[i]),nm,np)
  }
  }
  else
  {
    for(i in 1:length(auc)){
    np=length(which(labels==pos[i]))
    nm <- length(labels) - np
    pauc[i]<-pwilcox(nm*np*(1-auc[i]),nm,np)
    }
  }
  return(pauc)
}

pauclog <- function(auc,n=100,n.plus=0.5,labels=numeric(),pos=numeric()){
   pauc<-pauc(auc,n=n,n.plus=n.plus,labels,pos)
   pauclog<-array(0,length(auc))
   for(i in 1:length(auc)){
     pauclog[i]<-log10(pauc[i])
   }
  return(pauclog)
}


