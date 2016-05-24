defuzzify <-
function(XX){
 #function defuzzifies all elements of the list and returns vector of Steiner points
 #weighting measure is uniform one on [0,1], i.e. Lebesgue meassure
 temp_sum<-Msum(XX,0)
 if(is.null(temp_sum)==0){
 #subfun:
  defuzz<-function(X){
   #X is fuzzy set that will be defuzzified
   nl<-nrow(X)/2
   mids<-(X$x[1:nl]+X$x[(2*nl):(nl+1)])/2
   values<-mids[2:nl]+mids[1:(nl-1)]
   delta<-1/(nl-1)
   integral<-sum(values)*delta/2
   return(integral)
  }
  k<-length(XX)
    a<-rep(0,k)
    for (i in 1:k){
     a[i]<-defuzz(XX[[i]])
    }
  invisible(a)
  }
}
