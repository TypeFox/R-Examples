binaryTransform <- function(x, from, to){
  if (missing(from)){
    from <- sort(unique(c(unlist(x))))
  }
  if (length(from) != 2){
    stop("Binary data required")
  }
  
  x2 <- x
  x2[x==from[1]] <- to[1]
  x2[x==from[2]] <- to[2]
  x2
}

LinTransform <- function(graph, thresholds, from = c(0L, 1L), to = c(-1L, 1L), a, b)
{
  stopifnot(!missing(graph) & !missing(thresholds))

  if (missing(a) & missing(b))
  {  
    a <- (to[1]-to[2])/(from[1]-from[2])    
    b <- to[1] - a*from[1]
  }
  
  diag(graph) <- 0
  
  return(list(
    graph = graph/(a^2), 
    thresholds = thresholds/a - (b*rowSums(graph))/a^2))
}