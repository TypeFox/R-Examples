findNeighbour <- function(mom.mut,xspace,xpop,tested){

  new.pop <- NULL
  not.tested <- as.numeric(rownames(xspace))[! as.numeric(rownames(xspace)) %in% tested]

  for(i in 1:nrow(mom.mut)){
    term1 <- xspace[not.tested,]
    term2 <- matrix(as.matrix(mom.mut[i,]),nrow=nrow(term1),ncol=ncol(xspace),byrow=TRUE)
    dif <- apply((term1-term2)^2,1,sum)^(1/2)
    new <- as.numeric(names(which.min(dif)))
    new.pop <- c(new.pop,new)
    not.tested <- not.tested[! not.tested %in% new.pop]    
  }

return(new.pop)
}

