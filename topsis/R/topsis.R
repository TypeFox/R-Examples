topsis=function(decision = NULL, weights = NULL, impacts = NULL){
  if(missing(weights))
    stop("'weights' must be a numeric vector")
  if(missing(impacts))
    stop("'impacts' must be a character vector")
  if(! is.matrix(decision) |  is.data.frame(decision))
    stop("'decision' must be a matirx or data frame")
  if(length(weights) != ncol(decision))
    stop("length of 'weights' is not equal to number of columns")
  if(length(impacts) != ncol(decision))
    stop("length of 'impacts' is not equal to number of columns")
  if(! all(weights > 0))
    stop("weights must be positive numbers")
  if(! is.character(impacts))
    stop("impacts must be a character vector of '+' and '-' signs")
  if(! all(impacts == "+" | impacts == "-"))
    stop("impacts must be only '+' or '-' sign")
  weights <- weights/sum(weights)
  N <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for(i in 1:nrow(decision)){
    for(j in 1:ncol(decision)){
      N[i,j] <- decision[i,j] / sqrt(sum(decision[,j] ^ 2))
    }
  }
  W=diag(weights)
  V=N%*%W
  u <- as.integer(impacts == "+") * apply(V, 2, max) + 
    as.integer(impacts == "-") * apply(V, 2, min)
  l <- as.integer(impacts == "-") * apply(V, 2, max) + 
    as.integer(impacts == "+") * apply(V, 2, min)
  distance_u =function(x){
    sqrt(sum((x - u) ^ 2))
  }
  distance_l =function(x){
    sqrt(sum((x - l) ^ 2))
  }
  du <- apply(V, 1, distance_u)
  dl <- apply(V, 1, distance_l)
  score <- dl/(dl+du)
  return(data.frame(alt.row = 1:nrow(decision), score = score, rank = rank(-score)))
}
