
#Code for testing the probability of chance dimension restrictions

#the general idea is to take the box, remove a dimension restriction
#and see whether a box of size B' and size equal to the most restricted box could have been generated
#by chance given the mean inside the other bigger box.  

#this little bit of side code to evaluate box means when variables are removed:


nullprob <- function(dset, y=NULL, lbox){
  
  if(is.null(y)){   #If given null value for y, assume it is the last column of 
    y <- dset[,ncol(dset)] #the dataset.
  }
  
  if (length(lbox[[1]])>1){
    lvouts <- lvout(dset,y,lbox)
  } else{
    lvouts <- matrix(c(NA,NA,sum(y)/nrow(dset)),ncol=3)
    vect <- dset[,lbox[[1]][1]]
    invect <- (vect>lbox[[2]][1,1]) & (vect < lbox[[2]][1,2])
    attr(lvouts,"origtotin")  <- sum(invect) 
    attr(lvouts,"orighighin") <- sum(y[invect])
  }
  
  pvals <- vector(length=nrow(lvouts))
  
  for (i in 1:nrow(lvouts)){
  
    pbase <- lvouts[i,3]
    nbig  <- attr(lvouts,"origtotin")        #number of points in half restricted box
    intot <- attr(lvouts,"orighighin")       #number of high points in most restricted box
    #mean of box with restriction removed?
    
#Inverse framing:
    ptot <- 0
#    
#    for (j in 1:(intot-1)){
#    
#      ptot <- ptot + choose(nbig,j) * pbase^j * (1-pbase)^(nbig-j)
#      
#    }
#    
#    pvals[i] <- 1-ptot
   
#Test direct framing: (should come up with same results, except also work for intot=1) 

#Ok... this can just be calculated directly (and with better numerical 
#precision) using Rbeta, which uses pbeta -- 
#regularized incomplete beta function, which, with right parameters, functions
#as the upper cdf of the binomial dist.

#    for (j in nbig:intot){
#    
#      ptot <- ptot + choose(nbig,j) * pbase^j * (1-pbase)^(nchoosebig-j)
#      
#    }
#   
#    pvals[i] <- ptot
  #Rbeta(pbase,intot,nbig-intot+1)
   pvals[i] <- pbinom(intot-1, nbig, pbase, lower.tail=FALSE)
   
  }
  
  return(pvals)
  
}  



  
  
  
  
#    
#
#
#p1 = .2    
#p2 = .2
#
#pss <- c(1:20)
#
#for (k in 1:20){
#
#tot = 0
#
#n <- k
#
#for (i in 1:n){
#
#  ptot <- 0
#
#  for (j in 1:(i-1)){
#
#    ptot <- ptot + choose(n,j) * p2^j * (1-p2)^(n-j)
#
#  }
#  
#  tot <- tot + choose(n,i) * p1^i * (1-p1)^(n-i) * ptot
#  
#}
#
#pss[k] <- tot
#
#}
#
#
#for (i in 1:n){
#
#  ptot <- 0
#
#  for (j in 1:(i-1)){
#
#    ptot <- ptot + choose(n,j) * p2^j * (1-p2)^(n-j)
#
#  }
#  
#  tot <- tot + choose(n,i) * p1^i * (1-p1)^(n-i) * ptot
#  
#}
#
#
## I think this part is set up to test the lake model boxes
#
#n <- 74
#intot <- 60
#bigmean <- .7311
#
#ptot <- 0 
#
#for (i in 1:(intot-1)){
#
#  ptot <- ptot + choose(n,i) * bigmean^i * (1-bigmean)^(n-i)
#  
#}
#
#

