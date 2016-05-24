
#Code for evaluating box set when given in the same form as the output of boxext

#Read in something:
#setwd("C:/Documents and Settings/bryant/My Documents/algoOJT/ellipsoid")
#bar <- read.csv("2barrels-cross-3d.1000LH.csv", header=FALSE)



#penal <- function(n,hdisp=2,expo=2,vdisp=0){

#	return(pmax(0, (n-hdisp)^expo - vdisp))
  	
#}

#blist is constructed as follows:
#one element for each box
#blist[[i]][[1]] is a vector giving the dimension indices that are restricted
#blist[[i]][[2]] is a matrix with the restrictions on those dimensions


beval <- function(dset, blist, outdim="std", w = .5, dpparam=c(.05,2,2,0), bpparam=c(.01,1,1.5,0), justparams=FALSE){
  
  if(justparams){return(list(weta = w,dimpparams = dpparam, boxpparams = bpparam))}
  
  penal <- function(n,hdisp=2,expo=2,vdisp=0){

	 return(pmax(0, (n-hdisp)^expo - vdisp)) #WHY did I make this pmax?  3/22/07
  	
  }
  
  if (outdim=="std"){
    outdim <- ncol(dset)
  }   

  tpts <- nrow(dset)   #total points in the dataset
  this <- sum(dset[,outdim])  #total hi (onees) in the dataset
  
  dpscaler <- dpparam[1]
  bpscalar <- bpparam[1]
  
  B <- length(blist)
  
  bdims <- rep(0,B) #will be num of dimensions restricted for each box 
  
  for (i in 1:B){
    bdims[i] <- length(blist[[i]][[1]])
  }
  
  
# hisinset <- 0  #hi regret (one's) in the set
# tinset <- 0  #total point in the set
  
  #matrix of in/out vectors
  bincvecs <- matrix(TRUE,nrow=tpts,ncol=B)

  dnet <- c() #vector to keep track of all the dimensions restricted
  
  for (b in 1:B){
  
    dimvect <- blist[[b]][[1]]
    dnet <- c(dnet,dimvect)
    
    bmat <- blist[[b]][[2]]
  
    incvecs <- !logical(length=tpts)
  
    for (i in 1:length(dimvect)){
  		
      di <- dimvect[i]

      incvecs <-  incvecs & (dset[,di] >= bmat[i,1])
      incvecs <-  incvecs & (dset[,di] < bmat[i,2])
  
    } 
     
    bincvecs[,b] <- incvecs
    
    }
      

  
  masvecs <- logical(length=nrow(dset))
    
  for (b in 1:B){
    masvecs <- masvecs | bincvecs[,b]
  }

  dleft <- dset[masvecs,]

  hisinset <- sum(dleft[,outdim])
  tinset <- nrow(dleft)
  
#  fp <- (tinset-hisinset)/(tpts-this) #False positive rate
#  fn <- (this-hisinset)/(this) #False negative rate

#  dens <- 1-fp
#  cov  <- 1-fn
  
  dens <- hisinset/tinset
  cov  <- hisinset/this
   
  dimpenal <- dpscaler*sum(penal(bdims,dpparam[2],dpparam[3],dpparam[4]))
  
  bpenal <- bpscalar*penal(length(blist),bpparam[2],bpparam[3],bpparam[4])
  
  dct <- (dens)^w*(cov)^(1-w) #density/coverage tradeoff
  
  obj <- dimpenal + bpenal - dct
  
  bindim <- rep("-",outdim-1) 
  
  bindim[unique(dnet)] <- "X"  # Mark x's for dims that were restricted by the box set
  
  return(c(dens,cov,dct,obj,length(blist),min(bdims),max(bdims),sum(bdims),length(unique(dnet)),bindim))

}



