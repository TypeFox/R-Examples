components.graph <-
function(rl)
  {
	if (class(rl)!="landscape") 
  {
  stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
    cluster <- as.vector (rl$nodes.characteristics)["cluster"]
    b <- as.data.frame(table(cluster))
    c1 <- nrow(b)
    return(c1)
  }
