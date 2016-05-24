csigmap <-
function(n1,sd1,n2,sd2,weight=TRUE){
  #computes estimate of pooled sigma , see Cohen (1988, p. 66)
  # sd1 , sd 2 and n1, n2: standard deviation and sample size for group 1 and 2 respectively
      sigmap<-sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))      
   #cat("weighted sigmapooled equals " ,round(sigmap,dig=2))
      if (weight==FALSE){
   sigmap<-sqrt((sd1^2+sd2^2)/2)     #suggestion of Kraemer
    #cat("unweighted sigmapooled equals " ,round(sigmap,dig=2))
      }
  return(sigmap)}
