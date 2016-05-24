`olaps` <-
function(x,y,box.seq){

  nb <- length(box.seq)

  olap <- matrix(nrow=nb,ncol=nb)
  
  d <- ncol(x)
  n <- nrow(x)
  totint <- sum(y)
  
  for (i in 1:nb){
  
    pini <- in.box(x,box.seq[[i]]$box,d,boolean=TRUE)
  
    for (j in i:nb){
    
      pinj <- in.box(x,box.seq[[j]]$box,d,boolean=TRUE)
    
      supolap <- sum(pini & pinj)/n  #overlap of points
      covolap <- sum(y[pini & pinj])/totint  #overlap of interesting points
      
      olap[i,j] <- supolap
      olap[j,i] <- covolap
      
    }
  } 

  return(olap)
  
}

