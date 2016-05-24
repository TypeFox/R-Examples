adjustmeans <-
function(object, X, Y){
  
  newobject <- object
  justnodes <- newobject$nodes
  for (k in 1:length(justnodes)){
    ind <-  getsamples(justnodes[[k]],X,levelvec=justnodes$levelvec)
    indn <- (1:nrow(X))[- ind]
    attr(justnodes[[k]],"n") <- max(0.5,length(ind))
    attr(justnodes[[k]],"mean") <- if(length(ind)>0) mean(Y[ind])  else mean(Y)
  }
  newobject$nodes <- justnodes
  return(newobject)
}

