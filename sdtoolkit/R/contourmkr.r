
#This function takes a list containing matrices with coverage and density for 
#each box in a trajectory and translates it into a list of boxes each having 
#constant dimensions 


contourmkr <- function(margtraj){
  
  #marjtraj - as output by ranker()
  
  
  nboxes <- length(margtraj)  # How many boxes are in the sequence?
  
  bdims <- sapply(margtraj,nrow)  # How many dimensions are restricted in each box?
  
  maxd <- max(bdims) #What are the most dimensions restricted by any single box?
  
  transplist <- list()
  
  #for each box dimensionality
  for (d in 1:maxd){
  
    hasdim <- c(1:nboxes)[bdims >= d]   #which boxes have at least this many 
                                        #dimensions?
    
    transplist[[d]] <- matrix(ncol=2,nrow=nboxes)
    
    for (m in hasdim){ #for the boxes with this many dimensions or more
    
      transplist[[d]][m,] <- margtraj[[m]][d,3:2]   #extract the remove variable
                                                    #statistics for that number
                                                    #of dimensions (stored in  
                                                    #margtraj) 
    }
  
  }
  
  
#below is unnecessary, since plotting ignores the NA's anyway, and you can better
#use the identify function

#  for (i in 1:maxd){
#  
#    takeout <- apply(t(apply(transplist[[i]],1,is.na)),1,any)
#    
#    transplist[[i]] <- transplist[[i]][!takeout,]
#  
#  }
  
  return(transplist)

}
 