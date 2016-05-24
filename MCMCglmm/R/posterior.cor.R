"posterior.cor"<-function(x){

   if(is.null(dim(x))[1]){
     x<-as.matrix(x)
   }
   if(sqrt(dim(x)[2])%%1!=0){
    stop("cannot coerce rows of x into a matrix - the number of columns should be a triangular number")
   }
   coerce.cor<-function(x){cov2cor(matrix(x, sqrt(length(x)), sqrt(length(x))))}

   if(is.mcmc(x)==FALSE){
     warning("posterior.cor expecting mcmc object")
     t(apply(x, 1, coerce.cor))
   }else{   
     as.mcmc(t(apply(x, 1, coerce.cor)))
   }  
}

