systemGraphToGenerator <- function(g, failRate, repairRate) {
  numComp <- length(V(g))-2
  # Get igraph vertex IDs of non (s,t) nodes and (s,t) nodes
  C <- c()
  s <- 0
  t <- 0
  for(i in 1:(numComp+2)) {
    if(V(g)[i]$name == "s") {
      s <- as.numeric(V(g)[i])
      next
    }
    if(V(g)[i]$name == "t") {
      t <- as.numeric(V(g)[i])
      next
    }
    C <- c(C, as.numeric(V(g)[i]))
  }
  #print(C)
  #print(V(g)[C]$name)
  
  map <- list()
  for(state in (2^numComp-1):1) {
    bin <- digitsBase(state, ndigits=numComp)
    g2 <- delete.vertices(g, C[which(bin==0)])
    if(which(V(g2)$name=="t") %in% subcomponent(g2, which(V(g2)$name=="s"))) {
      #print(state)
      #print("works")
      map[[paste(bin, collapse="")]] <- list(failTo=c(), repairTo=c(), exitCount=0)
      cur <- paste(bin, collapse="") # current name
    } else {
      #print(state)
      #print("fails")
      next
    }
    for(i in 1:numComp) {
      if(bin[i]==0) {
        bin[i] <- 1
        map[[cur]]$repairTo <- c(map[[cur]]$repairTo, paste(bin, collapse=""))
        bin[i] <- 0
      } else {
        bin[i] <- 0
        g2 <- delete.vertices(g, C[which(bin==0)])
        if(which(V(g2)$name=="t") %in% subcomponent(g2, which(V(g2)$name=="s"))) {
          map[[cur]]$failTo <- c(map[[cur]]$failTo, paste(bin, collapse=""))
        } else {
          map[[cur]]$exitCount <- map[[cur]]$exitCount+1
        }
        bin[i] <- 1
      }
    }
  }
  
  Gnum <- matrix(0, nrow=length(map)+1, ncol=length(map)+1)
  Gchr <- matrix(0, nrow=length(map)+1, ncol=length(map)+1)
  const <- matrix(1, nrow=length(map)+1, ncol=length(map)+1)
  for(i in 1:length(map)) {
    for(k in map[[i]]$failTo) {
      j <- which(names(map)==k)
      Gnum[i,j] <- failRate
      Gchr[i,j] <- "F"
    }
    for(k in map[[i]]$repairTo) {
      j <- which(names(map)==k)
      Gnum[i,j] <- repairRate
      Gchr[i,j] <- "R"
    }
    if(map[[i]]$exitCount>0) {
      Gnum[i,length(map)+1] <- failRate*map[[i]]$exitCount
      Gchr[i,length(map)+1] <- "F"
      const[i,length(map)+1] <- map[[i]]$exitCount
    }
  }
  diag(Gnum) <- -rowSums(Gnum)
  list(G = Gnum, structure=list(G = Gchr, C = const))
}
