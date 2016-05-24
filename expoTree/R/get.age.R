# (recursive) get the age of a node in a phylogenetic 
# tree given the age of its ancestor
get.age <- function(tree,code,aage,tp=-Inf) {
  off <- which(tree$edge[,1] == code) 
  row <- which(tree$edge[,2] == code)
  root.len <- ifelse(is.null(tree$root.edge),0.0,tree$root.edge)
  edge.len <- ifelse(length(row)>0,tree$edge.length[row],root.len)
  parent   <- ifelse(length(row)>0,tree$edge[row,1],NA)
  age <- edge.len+aage
  if (length(off) > 0) {
    # internal node
    off.codes <- tree$edge[off,2]
    off.ages <- lapply(off.codes,get.age,tree=tree,aage=age,tp=tp)
    off.alive <- sapply(off.ages,function(x) any(x[,5]>0))
    oc <- do.call(rbind,off.ages)
    # check for extant tips
    #has.extant <- 1 %in% oc[,5]
    has.extant <- sum(off.alive)
    coord <- mean(sapply(off.ages,function(x) x[1,7]))
    ret <- cbind(parent,code,age,length(off),has.extant,aage,coord)
    return(rbind(ret,oc))
  } else {
    # tip
    coord <- NA
    if (! is.null(tree$coords)) {
      coord.row <- which(tree$coords[,1]==code)
      coord <- tree$coords[coord.row,2]
    }
    return(cbind(parent,code,age,length(off),age>tp,aage,coord))
  }
}

