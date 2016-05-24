# ' a summary-like function
# ' @param A coefficient vector
# ' @param labels name of the covariates
# ' @param X the dataset (named) if labels is null
# ' @param intercept boolean defining wether A contains an intercept or not
# ' @param ANOVA boolean to add Anova test for each coefficient
# ' @param print boolean to print ANOVA if computed
# ' @param Y the response variable if ANOVA is computed
# '@export
readY<-function(A=A,labels=NULL,X=NULL,intercept=TRUE,ANOVA=FALSE,print=FALSE,Y=NULL){
   if(is.null(labels)){
      labels=names(X)
   }
   if(intercept){labels=c("intercept",labels)}
   interp=cbind(A[A!=0],labels[A!=0])
   if(ANOVA & !is.null(X) & !is.null(Y)){
      Xred=X[,A[-intercept]!=0]
      reslm=lm(Y~.,data=data.frame(Xred))
      interp=summary(reslm) 
   }
   if(print){print(interp)}
   return(interp)
}