`predict.mdr` <- function(object, data=NULL, status=NULL, fold=NULL,...){
 
 TwoByTwo <- matrix(NA,ncol=2,nrow=2) 
 temp <- mdrEnsemble(object,data=data, new.status=status,fold=fold)
 res <- temp$result
 if(!is.null(status)) TwoByTwo <- table(res,status)
 output <- list(class=res,TwoByTwo=TwoByTwo,cv=temp$cv)
 class(output) <- "mdrPredict"
 output
} 
