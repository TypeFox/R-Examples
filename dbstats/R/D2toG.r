

  ########################
  #### function D2toG ####
  ########################


 D2toG <- function(D2,weights){
   if (class(D2)!="D2")
    stop(" 'D2' must be of class D2")
   
   n<-ncol(D2)
   
   if (missing(weights)||is.null(weights))
    weights<-rep(1,n)
   if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
   if (sum(weights<0)>0)
    stop("Weights array weights is not valid: some weights are negative")
   if (sum(weights)==0)
    stop("Weights array weights is not valid: sum(weights)=0")
  
  weights<-weights/sum(weights)
  
  
  G<-Gcalc(n,weights,D2)
  class(G)<-"Gram"
  return(G)

 }