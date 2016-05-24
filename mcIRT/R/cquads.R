### control supplied quads
cquads <- function(quads)
{

# NAMES  
  if(!all(as.vector(sapply(quads,function(e)names(e) == c("nodes","weights"))))){stop("\n quads: rename the supplied quads list!")}  

# NAs    
  anyNAs <- any(sapply(quads,function(e){
  any(any(is.na(e$nodes)),any(is.na(e$weights))  )
  }))
  
# length
  
  LEEEN <- lapply(quads,function(e){
    if(length(e$nodes) != length(e$weights)) {stop("quads: unequal lengths of weights and nodes!")}
    if(length(e$nodes) < 2){stop("quads: provide at least 2 nodes!")}
  })
  


  
# ZERO MEAN
  ergCTRL <- sapply(quads,function(x)
  {
  weighted.mean(x$nodes,x$weights)
  })
 
  if(!ergCTRL[1] < 0.001){stop("\n quads: Weighted sum of first node is not < 0.001. It is supposed to be zero!")}
  
  # VARIANCE = 1  
  
  ergCTRL2 <- sapply(quads,function(x)
  {
    sqrt(sum((x$nodes*x$weights - weighted.mean(x$nodes,x$weights))^2)/sum((x$nodes*x$weights)^2))
  })
  
  if(!ergCTRL2[1] > 0.99 & ergCTRL2[1] < 1.01){stop("\n quads: Weighted Standard Deviation of the first node is not 1!")}
  
  
  
}





