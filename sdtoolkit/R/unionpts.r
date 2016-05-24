

#I think this takes a dataset and a list of boxes, and finds all the points in
#the union of all the boxes.  


unionpts <- function(dset,blist){

  outdim <- ncol(dset)

  tpts <- nrow(dset)   #total points in the dataset
  this <- sum(dset[,outdim])  #total hi (onees) in the dataset


  B <- length(blist)

  bdims <- rep(0,B) #will be num of dimensions restricted for each box

  for (i in 1:B){
    bdims[i] <- length(blist[[i]][[1]])
  }

  #matrix of in/out vectors
  bincvecs <- matrix(TRUE,nrow=tpts,ncol=B)

  dnet <- c() #vector to keep track of all the dimensions restricted

  for (b in 1:B){

    dimvect <- blist[[b]][[1]]
    dnet <- c(dnet,dimvect)

    bmat <- blist[[b]][[2]]
    if(is.vector(bmat)){
      bmat <- matrix(bmat,ncol=2)
    }

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

  return(masvecs)
  
}
