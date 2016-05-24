mogavsToLinear <-
function(bestModel,y_ind,data,...){
  f2<-as.formula(paste(colnames(data)[y_ind],"~",paste(colnames(data)[which(bestModel==1)+1],collapse="+"),sep=""))
  modlm<-do.call("lm",list(formula=substitute(f2),data=substitute(data),model=TRUE,x=TRUE,y=TRUE))

  return(modlm)
}
