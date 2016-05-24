#' Calculate the weight matrix of a set of regimes on a phylogeny
#' 
#' These functions calculate weight matrices from regimes specified in phytools' simmap format.
#' \code{simmap.W} calculates the weight matrix for a set of regimes from a phylogeny
#' with a stored regime history. \code{.simmap.W} calculates the same matrix, but without checks and is 
#' generally run internally. 
#' 
#' @rdname simmap.W
#' @param tree either a tree of class "phylo" or a cache object produced by bayOU's internal 
#' functions. Must include list element 'maps' which is a simmap reconstruction of regime history.
#' @param pars a list of the parameters used to calculate the weight matrix. Only pars$alpha is
#' necessary to calculate the matrix, but others can be present.
#' 
#' @details \code{.simmap.W} is more computationally efficient within a mcmc and is used internally. The value
#' of \code{TotExp} is supplied to speed computation and reduce redundancy, and cache objects must be supplied as
#' the phylogeny, and the parameter \code{ntheta} must be present in the list \code{pars}.
#' @export
simmap.W <- function(tree,pars){
  if(class(tree)=="phylo"){
    X <- rep(NA,length(tree$tip.label))
    names(X) <- tree$tip.label
    cache <- .prepare.ou.univariate(tree,X)
  } else {cache <- tree}
  if(is.null(pars$ntheta)){
    pars$ntheta <- length(unique(names(unlist(cache$maps))))
  }
  nbranch <- length(cache$edge.length)
  maps <- cache$maps
  shifts <- unlist(lapply(maps,length),F,F)-1
  irow <- rep(1:nbranch,shifts+1)
  csbase <- cache$nH[irow]
  csmaps <- csbase+unlist(lapply(maps,cumsum),FALSE,TRUE)
  multips <- which(irow[2:length(irow)]==irow[1:(length(irow)-1)])
  csbase[multips+1] <- csmaps[multips]
  oW <- pars$alpha*(csbase-cache$height)
  nW <- (csmaps-csbase)*pars$alpha
  if(any(nW>500)){
    tmp <- ifelse(nW>500, exp(nW+oW), exp(oW)*(exp(nW)-1))
  } else {
    tmp <- exp(oW)*(exp(nW)-1)
  }
  bW <- matrix(0,nrow=nbranch,ncol=pars$ntheta)
  index <- irow + (as.integer(names(tmp))-1)*nbranch
  if(any(shifts>1)){
    tmp <- tapply(tmp,index,sum)
    bW[as.numeric(names(tmp))] <- tmp
  } else {bW[index] <- tmp}
  W=cache$branchtrace%*%bW
  W[,1] <- W[,1]+exp(-cache$height*pars$alpha)
  return(W)
}


.simmap.W <- function(cache,pars){
  nbranch <- length(cache$edge.length)
  maps <- cache$maps
  #Dangerous...may not have listed the shift (if shift occurs at node)
  shifts <- unlist(lapply(maps,length),F,F)-1
  #Index vector indicating which branch a given segment exists on
  irow <- rep(1:nbranch,shifts+1)
  #Height of the beginning of each branch
  csbase <- cache$nH[irow]
  #Height of the end of each segment
  csmaps <- csbase+unlist(lapply(maps,cumsum),FALSE,TRUE)
  #Determine which branches contain more than one segment
  multips <- which(irow[2:length(irow)]==irow[1:(length(irow)-1)])
  #Set segments with more than one shift per branch to start at end of last shift
  csbase[multips+1] <- csmaps[multips]
  #Exponential term 1
  oW <- pars$alpha*(csbase-cache$height)
  #Exponential term 2
  nW <- (csmaps-csbase)*pars$alpha
  #If value of expnential term is too large (resulting in overflow), then use approximation
  if(any(nW>500)){
    tmp <- ifelse(nW>500, exp(nW+oW), exp(oW)*(exp(nW)-1))
  } else {
    tmp <- exp(oW)*(exp(nW)-1)
  }
  #Set up branch weight matrix
  bW <- matrix(0,nrow=nbranch,ncol=pars$ntheta)
  #Set up index over matrix, so that values go to right row index and column, based on the name of the segment in maps
  index <- irow + (as.integer(names(tmp))-1)*nbranch
  if(any(shifts>1)){
    tmp <- tapply(tmp,index,sum)
    bW[as.numeric(names(tmp))] <- tmp
  } else {bW[index] <- tmp}
  W=cache$branchtrace%*%bW
  W[,1] <- W[,1]+exp(-cache$height*pars$alpha)
  return(W)
}

.parmap.W <- function(cache, pars){
  if(pars$k > 0){
    nbranch <- length(cache$edge.length)
    #create a vector called shifts that indicates the number of shifts for each branch
    nshifts <- table(pars$sb)
    shifts <- rep(0,nbranch)
    shifts[as.numeric(attributes(nshifts)$dimnames[[1]])]<- nshifts
    #Create an index equal to the number of segments that identifies the branch on which each segment is found
    irow <- rep(1:nbranch,shifts+1)
    #For now, starting height is just the height of the node
    csbase <- cache$nH[irow]
    #Calculate the ending height by sorting the edge.length and the location of shifts by their branch identity and location
    csadd <- c(cache$edge.length, pars$loc)
    tmp.o <- c(1:nbranch, pars$sb)
    names(csadd) <- tmp.o
    add.o <- order(tmp.o,csadd)
    csadd <- csadd[add.o]
    #Ending height of the segment
    csmaps <- csadd + csbase
    #We need to know what the ending theta is for each segment, so we sort pars$t2 as we did for pars$loc, but +1 because t2 is the ending regime
    t2index <- add.o[which(add.o > nbranch)]
    t2b <- c(rep(1,length(csmaps)))
    t2b[match(t2index,add.o)+1] <- pars$t2[t2index-nbranch]
    #Now we need to cascade these regime down the tree. We won't need to cascade sandwiches, as they are trapped on the branch they occur. So we find them below:
    loc.o <- order(pars$loc,decreasing=TRUE)
    sandwiches <- duplicated(pars$sb[loc.o])
    # And remove them:
    if(sum(sandwiches)>0){
      sb.down <- pars$sb[!sandwiches]
      t2.down <- pars$t2[!sandwiches]
    } else {sb.down <- pars$sb; t2.down <- pars$t2}
    #Now we order the sb's and t2's to prepare for a postorder tree traversal
    sb.o <- order(sb.down)
    sb.down <- sb.down[sb.o]
    t2.down <- t2.down[sb.o]
    sb.desc <- cache$bdesc[sb.down]
    #Loop traveling down the tree, saving all descendents that are from that shift into the vector censored. These branches cannot be modified by shifts further down the tree.
    censored <- NULL
    name.o <- names(csmaps)
    names(t2b) <- name.o
    for(i in 1:length(sb.desc)){
      sb.desc[[i]] <- sb.desc[[i]][!(sb.desc[[i]] %in% censored)]
      censored <- c(censored, sb.desc[[i]])
      t2b[name.o[name.o %in% sb.desc[[i]]]] <- t2.down[i]
    }
    names(csmaps) <- t2b
    multips <- which(irow[2:length(irow)]==irow[1:(length(irow)-1)])
    #Set segments with more than one shift per branch to start at end of last shift
    csbase[multips+1] <- csmaps[multips]
    #Exponential term 1
    oW <- pars$alpha*(csbase-cache$height)
    #Exponential term 2
    nW <- (csmaps-csbase)*pars$alpha
    #If value of expnential term is too large (resulting in overflow), then use approximation
    if(any(nW>500)){
      tmp <- ifelse(nW>500, exp(nW+oW), exp(oW)*(exp(nW)-1))
    } else {
      tmp <- exp(oW)*(exp(nW)-1)
    }
    #Set up branch weight matrix
    bW <- matrix(0,nrow=nbranch,ncol=pars$ntheta)
    #Set up index over matrix, so that values go to right row index and column, based on the name of the segment in maps
    index <- irow + (as.integer(names(tmp))-1)*nbranch
    if(any(duplicated(index))){
      tmp <- tapply(tmp,index,sum)
      bW[as.numeric(names(tmp))] <- tmp
    } else {bW[index] <- tmp}
    W=cache$branchtrace%*%bW
    W[,1] <- W[,1]+exp(-cache$height*pars$alpha)
  } else {
    W <- matrix(rep(1,cache$ntips),ncol=1)
  }
  return(W)
}

#' Calculate the weight matrix of a set of regimes on a phylogeny
#' 
#' These functions calculate weight matrices from regimes specified by a bayou formatted parameter list
#' \code{parmap.W} calculates the weight matrix for a set of regimes from a phylogeny
#' with a stored regime history. \code{.parmap.W} calculates the same matrix, but without checks and is 
#' generally run internally. 
#' 
#' @rdname parmap.W
#' @param tree either a tree of class "phylo" or a cache object produced by bayOU's internal 
#' functions. Must include list element 'maps' which is a simmap reconstruction of regime history.
#' @param pars a list of the parameters used to calculate the weight matrix. Only pars$alpha is
#' necessary to calculate the matrix, but others can be present.
#' 
#' @details \code{.parmap.W} is more computationally efficient within a mcmc and is used internally. 
#' @export
parmap.W <- function(tree, pars){
  if(class(tree)=="phylo"){
    X <- rep(NA,length(tree$tip.label))
    names(X) <- tree$tip.label
    cache <- .prepare.ou.univariate(tree,X)
  } else {cache <- tree}
  if(is.null(pars$ntheta)){
    pars$ntheta <- length(pars$theta)
  }
  if(pars$k > 0){
    nbranch <- length(cache$edge.length)
    #create a vector called shifts that indicates the number of shifts for each branch
    nshifts <- table(pars$sb)
    shifts <- rep(0,nbranch)
    shifts[as.numeric(attributes(nshifts)$dimnames[[1]])]<- nshifts
    #Create an index equal to the number of segments that identifies the branch on which each segment is found
    irow <- rep(1:nbranch,shifts+1)
    #For now, starting height is just the height of the node
    csbase <- cache$nH[irow]
    #Calculate the ending height by sorting the edge.length and the location of shifts by their branch identity and location
    csadd <- c(tree$edge.length, pars$loc)
    tmp.o <- c(1:nbranch, pars$sb)
    names(csadd) <- tmp.o
    add.o <- order(tmp.o,csadd)
    csadd <- csadd[add.o]
    #Ending height of the segment
    csmaps <- csadd + csbase
    #We need to know what the ending theta is for each segment, so we sort pars$t2 as we did for pars$loc, but +1 because t2 is the ending regime
    t2index <- add.o[which(add.o > nbranch)]
    t2b <- c(rep(1,length(csmaps)))
    t2b[match(t2index,add.o)+1] <- pars$t2[t2index-nbranch]
    #Now we need to cascade these regime down the tree. We won't need to cascade sandwiches, as they are trapped on the branch they occur. So we find them below:
    loc.o <- order(pars$loc,decreasing=TRUE)
    sandwiches <- duplicated(pars$sb[loc.o])
    # And remove them:
    if(sum(sandwiches)>0){
      sb.down <- pars$sb[!sandwiches]
      t2.down <- pars$t2[!sandwiches]
    } else {sb.down <- pars$sb; t2.down <- pars$t2}
    #Now we order the sb's and t2's to prepare for a postorder tree traversal
    sb.o <- order(sb.down)
    sb.down <- sb.down[sb.o]
    t2.down <- t2.down[sb.o]
    sb.desc <- cache$bdesc[sb.down]
    #Loop traveling down the tree, saving all descendents that are from that shift into the vector censored. These branches cannot be modified by shifts further down the tree.
    censored <- NULL
    name.o <- names(csmaps)
    names(t2b) <- name.o
    for(i in 1:length(sb.desc)){
      sb.desc[[i]] <- sb.desc[[i]][!(sb.desc[[i]] %in% censored)]
      censored <- c(censored, sb.desc[[i]])
      t2b[name.o[name.o %in% sb.desc[[i]]]] <- t2.down[i]
    }
    names(csmaps) <- t2b
    multips <- which(irow[2:length(irow)]==irow[1:(length(irow)-1)])
    #Set segments with more than one shift per branch to start at end of last shift
    csbase[multips+1] <- csmaps[multips]
    #Exponential term 1
    oW <- pars$alpha*(csbase-cache$height)
    #Exponential term 2
    nW <- (csmaps-csbase)*pars$alpha
    #If value of expnential term is too large (resulting in overflow), then use approximation
    if(any(nW>500)){
      tmp <- ifelse(nW>500, exp(nW+oW), exp(oW)*(exp(nW)-1))
    } else {
      tmp <- exp(oW)*(exp(nW)-1)
    }
    #Set up branch weight matrix
    bW <- matrix(0,nrow=nbranch,ncol=pars$ntheta)
    #Set up index over matrix, so that values go to right row index and column, based on the name of the segment in maps
    index <- irow + (as.integer(names(tmp))-1)*nbranch
    if(any(duplicated(index))){
      tmp <- tapply(tmp,index,sum)
      bW[as.numeric(names(tmp))] <- tmp
    } else {bW[index] <- tmp}
    W=cache$branchtrace%*%bW
    W[,1] <- W[,1]+exp(-cache$height*pars$alpha)
  } else {
    W <- matrix(rep(1,cache$ntips),ncol=1)
  }
  return(W)
}



