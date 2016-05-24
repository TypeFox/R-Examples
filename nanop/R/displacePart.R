displacePart <- function(nanop, sigma=NA,
                         rcenter=FALSE,
                         latticep=4.08) {
  onanop <- nanop 
  if(length(nanop) < 1)
    return() 
  if(rcenter)
    if(is.null(attributes(nanop)$center))
      center <- runif(3, max=max(latticep))
    else
      center <- attributes(nanop)$center
	  
  nAtomTypes <- attributes(nanop)$nAtomTypes
  
  if (is.na(sigma[1]))
    sigma <- attributes(nanop)$sigma
  if (is.na(sigma[1]))
    sigma <- 0
  if (length(sigma) < nAtomTypes){
    warning("vector <sigma> length is less than number of atom types nAtomTypes\n", immediate. = TRUE)
	sigma <- rep(sigma[1], nAtomTypes)
  }
  
  atomType <- attributes(nanop)$atomType
  tt <- 1:max(atomType)
  if(min(atomType) < 0)
    tt <- c(tt, -1:min(atomType))
  nAt <- 0 # number of atoms of each type
  for(i in 1:nAtomTypes){
    nAt[i] <- length( which(atomType == tt[i]) )  
  }
  
  iAtomType <- list() # positions of atoms in whole array 'nanop' for each atom type
  for(i in 1:nAtomTypes){
    iAtomType[[i]] <-  which(atomType == tt[i])
  }
  
  good <- 1 #have all sigma needed
  for(i in 1:nAtomTypes){
    if(is.na(sigma[i]))
	  good <- 0
  }
  
  if( nAtomTypes == 1 )  {
    ## uniform
    np <- rnorm(nrow(nanop)*3, mean = 0, sd = sqrt(sigma[1]))
    nanop <- nanop + np
    atomType_new <- atomType
  }
  else if( good==1 && nAtomTypes >1 ) {
    nanop_new <- matrix(nrow=0, ncol=3)
    atomType_new <- 0
    for(i in 1:nAtomTypes){
	    if( nAt[i] >= 1 ){
        nanop_tmp <- nanop[iAtomType[[i]], ]  	  
	      if(nAt[i]==1) ## one atom of given type
          nanop_tmp <- as.matrix(t(nanop_tmp))
		   
	      np <- rnorm(nrow(nanop_tmp)*3, mean = 0, sd=sqrt(sigma[i])) 
        nanop_tmp <- nanop_tmp + np
		    nanop_new <- rbind(nanop_new, nanop_tmp)
        atomType_new <- c(atomType_new, rep(tt[i], nAt[i]))
	    }
	  }
    nanop <- nanop_new
    atomType_new <- atomType_new[-1]
  }
  else 
    stop("Non-sensible input argument")
    
  attributes(nanop) <- attributes(onanop)
  attributes(nanop)$atomType <- atomType_new
  nanop
}
