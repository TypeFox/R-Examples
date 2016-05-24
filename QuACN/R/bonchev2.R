bonchev2 <- function(g, dist=NULL, wien=NULL){
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  if(is.null(wien)){
    wien <- wiener(g)
  }
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }

  rho <- max(dist)
  ki <- table(dist)[2:(rho+1)]
  i <- as.numeric(names(ki))
  In <- i*ki * log2(i)

  wien*log2(wien)-sum(In)
}
