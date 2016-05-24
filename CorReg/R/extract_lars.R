# ' extrait un des mod?les du lars
extract_lars<-function(lars=lars,qui=qui,X=X,Y=Y){
  coef=coef(lars)[qui,]
  qui=which(coef!=0)
  formule_opt=paste("Y~",paste(names(qui),collapse="+"))
  base_loc=data.frame(X[,qui+1],Y)
  names(base_loc)=c(names(qui),"Y")
  mod_opt=lm(formule_opt,data=base_loc)
  return(mod_opt)
}