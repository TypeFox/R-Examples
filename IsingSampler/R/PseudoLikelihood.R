# Pseudo likelihood:
IsingPL <- function(
  x, # Vector or data frame containing data
 graph, thresholds, beta, responses = c(0L,1L))
{
  stopifnot(isSymmetric(graph))  
  stopifnot(length(responses)==2)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  N <- nrow(graph)
  
  # If x is vector, turn it into a matrix:
  if (is.vector(x)){
    x <- t(x)
  }
  
  if (is.data.frame(x)){
    x <- as.matrix(x)
  }
  
  stopifnot(is.matrix(x))
  
  # Compute pseudo likelihood:
  PL <- PseudoLikelihood(x, graph, thresholds, beta, responses, logis = TRUE)
  
  
  return(PL)
}

# Optimisation:
EstimateIsingPL <- function(data, responses, beta = 1, ...){
  if (missing(responses)){
    responses <- sort(unique(c(data)))
  }

  if (length(responses) != 2){
    stop("Binary data required")
  }
  
  # Optimisation function:
  optimFun <- function(par, data){
    Ni <- ncol(data)
    stopifnot(length(par) == (Ni*(Ni-1)/2)+Ni)
    graph <- matrix(0,Ni, Ni)
    graph[upper.tri(graph)] <- par[1:(Ni*(Ni-1)/2)]
    graph[lower.tri(graph)] <- t(graph)[lower.tri(graph)]
    thresholds <- par[(Ni*(Ni-1)/2+1):length(par)]
    -2*IsingPL(data, graph, thresholds, 1, responses = sort(unique(c(data))))
  }
  
  # Run optimizer:
  Ni <- ncol(data)
#   optimRes <- optim(rep(0,(Ni*(Ni-1)/2)+Ni), optimFun, data = Data, ...)
  optimRes <- nlm(optimFun, rep(0,(Ni*(Ni-1)/2)+Ni), data = data, ...)
  
  # Cunstruct graph and thresholds:
  par <- optimRes$estimate
  graph <- matrix(0,Ni, Ni)
  graph[upper.tri(graph)] <- par[1:(Ni*(Ni-1)/2)]
  graph[lower.tri(graph)] <- t(graph)[lower.tri(graph)]
  thresholds <- par[(Ni*(Ni-1)/2+1):length(par)]
  
  return(list(
    graph = graph,
    thresholds = thresholds,
    results = optimRes))
}