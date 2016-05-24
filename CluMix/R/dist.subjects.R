dist.subjects <-
function(data){
#function(data, type=list()){
# !! to be done: allow also asymmetric binary variables
  
  # if all variables are numeric, use Euclidean distance
  dc <- sapply(data, data.class)
  if(all(dc == "numeric"))
    D <- dist(data)
  
  # if not, use Gower's distance with Podani's extension
  else{
    # !! depending on type, define asymmetric binary variables for parameter asym.bin  
  
    # binary variables have to be numeric
    K <- sapply(data[,dc == "factor", drop=FALSE], function(x) length(levels(x)))
    bin <- names(K)[K == 2]
    data[,bin] <- sapply(data[,bin], function(x) as.numeric(x) - 1)
  
    D <- FD::gowdis(x=data, ord="metric") # asym.bin=!!
  }
  return(D)
}
