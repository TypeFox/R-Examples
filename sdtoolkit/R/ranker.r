`ranker` <-
function(checpts, x, y,infolist,npts,ninter){

#checpts: a vector of which boxes in the trajectory to evaluate
#x:  The current (covered) input set
#y:  The current (covered) output vector
#infolist:  as returned by  traj.info - list with box elements and box stats

#RETURNS: A list of length checpts, each component of which is a matrix
#each matrix has the column number [CHECK!], density, coverage and support of 
#remove variable statistics for each box

  d <- ncol(x)
  n <- nrow(x)
  margsum <- sum(y)
  
  eibbox <- list(length=length(checpts)) #Extra Info By Box 
  #print(cat(checpts,"bah"))
  #flush.console()
  for(h in 1:length(checpts)){
  
    ### CODE REMNANT 1 WAS HERE, MOVED TO END OF FILE JUST IN CASE


    box.curr <- infolist$box[[checpts[h]]]   #get current box definition

    rdims <- c(1:d)[infolist$dimlist[[checpts[h]]]$either] #which inds are
    currentdimthing <- rdims

    ### CODE REMNANT 2 WAS HERE

    stats=matrix(nrow=length(rdims), ncol=4)

    ranki <- 0
    
    while(length(currentdimthing)>1){
    
      ranki <- ranki+1       
      
      supstore <- meanstore <- covstore <- istore <- vector(length=length(currentdimthing)-1)
      
      for (i in currentdimthing){
        
        trdims <- currentdimthing[currentdimthing!=i]  #restrict all dimensions except the current one, to assess the 'remove variable statistics'
        
        x.ind.curr <- rep(TRUE,n)
          
        for (j in trdims){
        
          x.ind.curr <- x.ind.curr & (x[,j]>= box.curr[1,j]) & (x[,j] <= box.curr[2,j])
    
          x.curr <- x[x.ind.curr,] # & x.ind,]  #why have the  &x.ind??

    
#          xy.list$x[[k]] <- x.curr
       
        }
 
        
        box.mass.curr <- sum(x.ind.curr)/npts
 
        ssi <- c(1:length(currentdimthing))[currentdimthing==i]
 
        istore[ssi] <- i
        supstore[ssi]  <- box.mass.curr
        meanstore[ssi] <- mean(y[x.ind.curr])
        covstore[ssi]  <- sum(x.ind.curr)*meanstore[ssi]/ninter  
 #       covstore[ssi]  <- supstore[ssi]*meanstore[ssi]*n/ninter
        
      }
      
      #identify which i in currentdimthing reduces the mean the least
      
      maxi <- which.max(meanstore)
      
      leastimp <- currentdimthing[maxi]
 
      #store that eye and the stats somewhere, remove that i from current dimthing
      
      stats[ranki,] <- c(leastimp,meanstore[maxi],covstore[maxi],supstore[maxi]) 
      
#      extralist$ranking[ranki] <- leastimp #currentdimthing[istore[leastimp]]
      
      currentdimthing <- currentdimthing[-maxi]
      
    }
    
    #cat(dim(stats),ranki)  #Added 2009-10-25 as diagnostic
    stats <- matrix(stats[(ranki+1):1,],ncol=4)
    
    gmean  <- margsum/n
    marcov <- margsum/ninter          #will need to modify for boxes farther in coverage seq
    gsup   <- n/npts           #will need to modify also
      
    stats[1,] <- c(currentdimthing,gmean,marcov,gsup)
      
    eibbox[[h]] <-  stats
    
  }
  
  return(eibbox)
  
}

### CODE REMNANT 1 ####

#    box.curr <- box.seq$box[[checpts[h]]]

#     
#  #only bother with dims that have been restricted
#    dimr <- dimchecker(x,box.curr)
#    rdims <- c(1:d)[dimr$either]
#    currentdimthing <- rdims
#
#   identical(infolist$box[[h]],box.curr)


#### CODE REMNANT 2:
   
#  extralist$stats[[h]] <- vector(length=length(rdims))
# r\turn(extralist)
#}
#  x.ind <- rep(TRUE, n)
#  xy.list <- list()
#  