best.layout <-
function(x){
  
  x1 <- merge(1:x,1:x) # until a graph with 9,000 spaces
  x2 <- abs(apply(x1,1,function(x,des){x <- unlist(x);y <- (x[[1]]*x[[2]])-des;return(y)},des=x))
  vw <- x1[which(x2 <2),] # which acoommodation gives the least empty space in graph
  x3 <- abs(apply(vw,1,function(x){y <- abs(x[1]-x[2]); return(y)}))
  x4 <- vw[which(x3 == min(x3))[1],] # best dimensions for the layout
  res <- as.vector(unlist(x4))
  return(res)
}
