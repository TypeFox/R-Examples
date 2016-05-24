"posterior.evals"<-function(x){

   if(is.null(dim(x))[1]){
     x<-as.matrix(x)
   }

   if(sqrt(dim(x)[2])%%1!=0){
    stop("cannot coerce rows of x into a matrix - the number of columns should be a triangular number")
   }

   coerce.eval<-function(x){eigen(matrix(x, sqrt(length(x)), sqrt(length(x))))$values}

   if(is.mcmc(x)==FALSE){
     warning("posterior.evals expecting mcmc object")
     t(apply(x, 1, coerce.eval))
   }else{   
     as.mcmc(t(apply(x, 1, coerce.eval)))
   }  
}


