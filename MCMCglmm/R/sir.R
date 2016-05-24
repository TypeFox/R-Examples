"sir"<-function(formula1=NULL, formula2=NULL){

   if(class(formula1)!="formula"){stop("formula not passed to formula1 in sir")}
   if(class(formula2)!="formula"){stop("formula not passed to formula2 in sir")}

   formula1<-update.formula(formula1, ~.-1)
   formula2<-update.formula(formula2, ~.-1)

 #  X1<-sparse.model.matrix(formula1)
 #  X2<-sparse.model.matrix(formula2)

   X1<-model.matrix(formula1)
   X2<-model.matrix(formula2)

   if(dim(X1)[2]!=dim(X2)[2]){
     stop("sir formulae invalid: factor levels of intersecting variables have to be the same")
   }else{
     X1%*%t(X2)
   }
}


