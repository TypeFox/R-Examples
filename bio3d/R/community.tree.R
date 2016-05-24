community.tree <- function(x, rescale=FALSE){

  ## Check for presence of igraph package
  oops <- requireNamespace("igraph", quietly = TRUE)
  if (!oops) {
     stop("igraph package missing: Please install, see: ?install.packages")
  }

  if(class(x) != "cna"){
    stop("Input should be a 'cna' class object as obtained from cna()")
  }

  rescaling <- function(membership){
    original.comms <- unique(membership)
    new.comms <- c(1:length(unique(membership)))
      
    a<-1  ## index to keep track of community number
    for(j in 1:length(membership)){
      membership[membership == original.comms[a]] <- new.comms[a]
      a <- a +1
    }
    return(membership)
  }
  
  num.of.nodes <- length(igraph::V(x$network))
  membership <- c(1:num.of.nodes)

  merge.table <- x$communities$merges

  membership.table <- matrix(NA, nrow=dim(merge.table)[1]+1, ncol=num.of.nodes)

  membership.table[1,] <- c(1:num.of.nodes)
  num.of.comms <- c(num.of.nodes, rep(NA,dim(merge.table)[1]))
  
  for(i in 1:dim(merge.table)[1]){
    comm.number <- num.of.nodes + i

    if(merge.table[i,1] < num.of.nodes){
      membership[merge.table[i,1]] <- comm.number
    }
    else{
      change.inds <- which(membership == merge.table[i,1])
      membership[change.inds] <- comm.number
    }
 
    if(merge.table[i,2] < num.of.nodes){
      membership[merge.table[i,2]] <- comm.number
    }
    else{
      change.inds <- which(membership == merge.table[i,2])
      membership[change.inds] <- comm.number
    }

    membership.table[i+1,] <- membership  ## i+1 because the first line is where each node forms a separated community (it will match also the modularity values)
    num.of.comms[i+1] <- length(unique(membership))
  }

  ## Rescale community number starting from 1
  if(rescale){
    membership.table <- t(apply(membership.table,1,rescaling))
  }
  
  output <- list("tree" = membership.table,
                 "modularity" = x$communities$modularity,
                 "num.of.comms" = num.of.comms)

  return(output)
}
