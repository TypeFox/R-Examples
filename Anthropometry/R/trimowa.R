trimowa <- function(data,w,numClust,alpha,niter,algSteps,ah=c(23, 28, 20, 25, 25),verbose){

 asig <- c() ; copt <- c() ; s <- c(0) ; n <- c() ; no.trim <- c()       
 res_qq <- c() ; difs <- c() ; bh <- c() ; bl <- c() ; trim <- c()

  n1 <- nrow(data)
  n2 <- floor(n1 * (1 - alpha))
  s = s + (n1 - n2) #number of trimmed individuals.
 
  #Constants that appear in "An Optimisation Approach to Apparel Sizing" of McCulloch et al:
  bh <- (apply(as.matrix(log(data)), 2, range)[2,] - apply(as.matrix(log(data)), 2, range)[1,]) / ((numClust-1) * 8) 
  bl <- -3 * bh
  ah <- ah
  al <- ah / 3

  #Data processing:
  num.persons <- dim(data)[1]
  num.variables <- dim(data)[2]
  datam <- as.matrix(data)
  datat <- aperm(datam, c(2,1))                     
  dim(datat) <- c(1, num.persons * num.variables)   
  rm(datam)

  #Computing the dissimilarity matrix:
  D <- getDistMatrix(datat, num.persons, num.variables, w, bl, bh, al, ah, verbose)
  rm(datat)

  n <- dim(D)[1] #number of individuals in each bust class.
  no.trim <- floor(n * (1 - alpha)) #number that left after trimming.

  #"Garbage collector" to release as much memory as possible:
  gc(verbose = FALSE)

   #Trimmed K-medoids:
   results <- trimmedoid(D, numClust, alpha, niter, algSteps, verbose)

    #asig gives the cluster to which each individual belongs and copt are the centroids of the clusters:
    asig <- results$asig
    #Obtain the current medoids obtained for each bust class, regarding the whole database:
    copt <- as.numeric(rownames(data)[results$copt])

    #qq is the vector that only contains the individuals after the trimmed procedure. 
    #With the following commands, we obtain the current trimmed individuals in the whole database:
    res_qq <- results$qq 
    difs <- setdiff(1:n,res_qq) 
    trim <- as.numeric(rownames(data[difs, ]))

 return(list(cases=copt,numTrim=s,numClass=n,noTrim=no.trim,C1=bh,C2=bl,C3=ah,C4=al,asig=asig,discarded=trim))
}
