Bvar <-
function(XX,theta=1/3){
  #calculates the variance of k polygonal fuzzy numbers with same levels
  #if necessary just use translator first to assure same alpha levels
  #theta ... is weight in the def of the bertoluzza metric
  sample_mean<-Mmean(XX,0)
  if(is.null(sample_mean)==0){
   if(length(XX)==1){
    v<-0
    }
  if(length(XX)>=2){
   k<-length(XX)
     temp<-rep(0,k)
     for (i in 1:k){
      temp[i]<-(bertoluzza(XX[[i]],sample_mean,theta,0))^2
      }
    v<-mean(temp)
   }
  invisible(v)
 }
}
