Media_pvalue<-function(iTest,jTest,times,selectvec){

  pout<-permute(iTest,jTest,times,selectvec)
  lengthtemp<-length(selectvec)
  result<-c(1:lengthtemp)
  for(i in 1:length(selectvec)){ result[i]<-pnorm(dis_value(iTest,jTest,selectvec[i]),pout[[1]][[i]],pout[[2]][[i]])
  }
  return(result)
}
