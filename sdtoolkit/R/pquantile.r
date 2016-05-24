
#quantile finding function that deals with repeated values in a way that 
#is appropriate for PRIM as described below

pquantile <- function(x,alpha,ptype){

#for PRIM, if the quantile is in between identical values, need to narrow it
#means raising threshold for low end, lowering for high end.
#ptype indicates which we are lowering or raising. 
#Eg, if we are just counting entries, the value associated with the first 
#quartile of this set of numbers is 2:  c(1,2,2,4,6,10,11,13). 
#However, we do want to make sure we cut off the first 3 entries, not just the 
#first two... [though... why wouldn't entry>2 take care of that... --
#oh, would only be issue if using >= instead of >?  Which... is not the case in
#peel.one, so... what?

#alpha patience/quantile parameter
  
  N <- length(x)          
  
#  if(N==45){browser()}                 
 # print(N) #added as diagnostic 2009-10-25
  xs <- sort(x)
  
  intm <- alpha*N
 # cat(N,alpha,intm,floor(intm),ceiling(intm))
  
  #pthresh <- .5*(xs[floor(intm)]+xs[ceiling(intm)])#Commented out 2009-10-25
  #provisional threshold, assuming all x values are unique:
  pthresh <- .5*(xs[max(floor(intm),1)]+xs[max(ceiling(intm),1)])
  #The "max" function is for when using really small peel.alpha's, which 
  #can lead to intm < 1.
  
  #check for lots o discrete values
 # print(N)
 # print(intm)
 # browser()
 #Deal with special cases:
  if (xs[max(floor(intm),1)]==xs[max(ceiling(intm),1)] & 
			min(xs)!=max(xs)){ #when min(xs)==max(xs), the warnings will be triggered
			                   #but the candidate threshold should already be valid
		
    if ((floor(intm) <= 1 | ceiling(intm) >= N)){
    
      if (ptype=="lowend"){ #keep raising index up
    
        #if (xs[floor(intm)] == xs[floor(intm)-1]){ 
          # means have repeated values at the quantile messing up pasting until number is different
           
          i <- 2
          
          while ((i<=N) & (xs[i]==xs[i-1])){
            i <- i+1
          }
          
          if(i <= N){
            pthresh <- .5*(xs[i] + xs[i-1])
          } else {
						pthresh <- NA
					}
        
      } 
      
      if (ptype=="highend"){# keep lowering index until number is different
        
        #if (xs[floor(intm)] == xs[floor(intm)+1]){ 
          # means have repeated values at the quantile messing up pasting until number is different
        
          i <- N-1
          
          while ((i > 1) & (xs[i]==xs[i+1])){
            i <- i-1
          }
          
          if (i > 0){
            pthresh <- .5*(xs[i] + xs[i+1])
          } else {pthresh <- NA}
       } 
      } else if (ptype=="lowend"){ #keep raising index up
    
        if (xs[floor(intm)] == xs[floor(intm)-1]){ 
          # means have repeated values at the quantile messing up pasting until number is different
           
          i <- ceiling(intm)
          
          while ((i<=N) & (xs[i]==xs[i-1])){
            i <- i+1
          }
          
          if(i <= N){
            pthresh <- .5*(xs[i] + xs[i-1])
          } else {pthresh <- NA}
        }
      } else if (ptype=="highend"){# keep lowering index until number is different
        
        if (xs[floor(intm)] == xs[floor(intm)+1]){ 
          # means have repeated values at the quantile messing up pasting until number is different
        
          i <- floor(intm)
          
          while ((i > 1) & (xs[i]==xs[i+1])){
            i <- i-1
          }
          
          if (i > 0){
            pthresh <- .5*(xs[i] + xs[i+1])
          } else {pthresh <- NA}
        }
    }
  }

  return(pthresh)

}











