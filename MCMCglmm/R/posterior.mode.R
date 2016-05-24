"posterior.mode"<-function(x, adjust=0.1, ...){

   if(is.mcmc(x)==FALSE){
     warning("posterior.mode expecting mcmc object")
   }

  find.mode<-function(x,adjust,...){
    dx<-density(x, adjust=adjust, ...)
    dx$x[which.max(dx$y)]
  }
  apply(as.matrix(x), 2, find.mode, adjust=adjust, ...)
}

