Ddivergence<-function(CA=NULL, CB=NULL, n=10000){

  if(requireNamespace("mvtnorm", quietly = TRUE)==FALSE){stop("mvtnorm not loaded")}

  if(dim(CA)[1]!=dim(CB)[1] | dim(CA)[2]!=dim(CB)[2] | dim(CA)[1]!=dim(CA)[2]){
     stop("matrices must be the same dimension and square")
  }

  xi<-mvtnorm::rmvnorm(n, rep(0,dim(CA)[1]), CA)
  fx<-mvtnorm::dmvnorm(xi, rep(0,dim(CA)[1]), CA)
  gx<-mvtnorm::dmvnorm(xi, rep(0,dim(CA)[1]), CB)

  mean(sqrt(0.5*((fx-gx)^2)/(fx+gx)))

}

