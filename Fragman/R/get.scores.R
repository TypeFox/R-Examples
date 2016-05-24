get.scores <-
function(my.scores, mark="mark"){
  y <- lapply(my.scores, function(x){x$wei})
  n <- lapply(my.scores, function(x){length(x$wei)})
  nn <- max(unlist(n))
  da <- data.frame(matrix(NA, byrow=T, ncol=nn, nrow=length(my.scores)))
  rownames(da) <- names(my.scores)
  for(i in 1:dim(da)[1]){
    da[i, 1:(length(y[[i]]))] <- y[[i]]
  }
  col.na <- vector(mode="character")
  for(k in 1:nn){
    col.na[k] <- paste(mark,"A.",k,sep="")
  }
  names(da) <- col.na
  return(da)
}
