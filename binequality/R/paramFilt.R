paramFilt <-
function(paramsList, fitComb){
  for(i in names(paramsList)){
    undef.i<-numeric(0)
    params.i<-paramsList[[i]]
    if(i != 'GG'){
      if(i == 'DAGUM'){
        params.i$tau<-1
      }#end if 
      if(i == 'BETA2'){
        params.i$sigma<-1
      }#end if 
      if(i == 'SINGMAD'){
        params.i$nu<-1
      }#end if 
      if(i == 'LOGLOG'){
        params.i$tau<-1
        params.i$nu<-1
      }#end if 
      if(i == 'PARETO2'){
        params.i$tau<-1
        params.i$nu<-1
      }#end if 
      undef.mu<-which(-params.i$nu>1/params.i$sigma|1/params.i$sigma>params.i$tau)
      undef.var<-which(-params.i$nu>2/params.i$sigma|2/params.i$sigma>params.i$tau)
      undef.i<-c(undef.mu,undef.var)
      undef.i<-unique(undef.i)
    }#end if not GG
  if(length(undef.i)>0){
    useDist.i<-which(fitComb$distribution==i)
    fitComb[useDist.i,c('estMean','var','cv','cv_sqr','gini','theil','median','sd')][undef.i,]<-NA
  }#end if length(def.i)
}#end for i
  return(fitComb)
}
