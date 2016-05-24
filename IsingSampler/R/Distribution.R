IsingSumLikelihood <- function(graph, thresholds, beta, responses = c(0L,1L))
{
  stopifnot(isSymmetric(graph))  
  stopifnot(length(responses)==2)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  N <- nrow(graph)
  Allstates <- do.call(expand.grid,lapply(1:N,function(x)c(responses[1],responses[2])))
  P <- exp(- beta * apply(Allstates,1,function(s)H(graph,s,thresholds)))
  SumScores <- rowSums(1*(Allstates==1))

  df <- ddply(data.frame(Sum = SumScores, P = P),"Sum",summarize,P=sum(P))
  df$P <- df$P / sum(df$P)
  return(df)
}

IsingLikelihood <- function(graph, thresholds, beta, responses = c(0L,1L))
{
  stopifnot(isSymmetric(graph))  
  stopifnot(length(responses)==2)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  N <- nrow(graph)
  Allstates <- do.call(expand.grid,lapply(1:N,function(x)c(responses[1],responses[2])))
  P <- exp(- beta * apply(Allstates,1,function(s)H(graph,s,thresholds)))
  df <- cbind(Probability = P / sum(P), Allstates)
  return(df)
}

IsingStateProb <- function(s,graph,thresholds,beta,responses=c(0L,1L))
{
  if (!is.list(s)) s <- list(s)
  N <- length(s[[1]])
  Allstates <- do.call(expand.grid,lapply(1:N,function(x)responses))
  Dist <- exp(- beta * apply(Allstates,1,function(s)H(graph,s,thresholds)))  
  Z <- sum(Dist)  
  
  sapply(s, function(x)exp(-beta*H(graph,x,thresholds))/Z)
}