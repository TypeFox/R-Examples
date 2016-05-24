get_endnodes <-
function(info){
  
  endnodes      <- list()
  endnodes[[1]] <- 1
  for(j in 1:nrow(info)){
    endnodes[[j+1]] <- numeric(length=(j+1))
    what <- c(info[j,"left"],info[j,"right"])
    delete <- info[j,"number"]
    where  <- which(endnodes[[j]]==delete)
    endnodes[[j+1]][c(where,where+1)] <- what
    endnodes[[j+1]][-c(where,where+1)] <- endnodes[[j]][-where]
  }
  endnodes <- endnodes[[length(endnodes)]]
  return(endnodes)
  
}
