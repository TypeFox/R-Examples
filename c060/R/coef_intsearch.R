
###########################################################################################################
coef.sum.intsearch<-function(object,...){
  # get coef for a object from fit object after running interval search 
   
  f1 <- object$cvreg
  res<-f1$glmnet.fit
  cof <- as.vector(coef(res, s=object$lambda))
  names.cof <- rownames(res$beta)
  cofn <- cof[which(cof != 0)]
  names(cofn) <- names.cof[which(cof != 0)] 
  bet <- res$beta[match(names(cofn), rownames(res$beta)),]
  
  return(cofn) 
}  

