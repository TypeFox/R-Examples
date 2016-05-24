ei.reg.bayes <- function(formula, data, sample=1000, weights 
= NULL, truncate=FALSE) { 
  D <- model.frame(formula, data = data)
  G <- D[[2]]
  T <- D[[1]]

  idx.r <- apply(G,1,sum)
  idx.c <- apply(T,1,sum)

  countG <- countT <- propG <- propT <- FALSE

  if(all(as.integer(T)==T) && all(T)>=0){
    countT <- TRUE
  }
  else{
    propT <- TRUE
  }
  
  if(all(as.integer(G)==G) & all(G)>=0){
    countG <- TRUE
  }
  else{propG <- TRUE}
  
  if(countT & countG){
    if(all(idx.r == idx.c)){
      G <- G/idx.r
      T <- T/idx.c
      countT <- countG <- FALSE
      propT <- propG <- TRUE
    }
    else{
      stop("row and column count totals unequal in some precincts -
please respecify data")
    }
  }
  
  if(countT & propG){
    if(!all(0 <= G && G <= 1)){
      stop("row proportions are not within [0,1] - please respecify
data")
    }
    else{
      T <- T/idx.c
      propT <- TRUE
      countT <- FALSE
    }
  }

  if(propT & countG){
    if(!all(0 <= T && T <= 1)){
      stop("column proportions are not within [0,1] - please respecify
data")
    }
    else{
      G <- G/idx.r
      propG <- TRUE
      countG <- FALSE
    }
  }

  if(propT & propG){
    idx.r <- apply(G,1,sum)
    idx.c <- apply(T,1,sum)

    if(!all(round(idx.r, digits=3)==1)){
      stop("row marginals are proportions that do not sum to 1 - please
respecify data")
    }

    if(!all(round(idx.c, digits=3)==1)){
      stop("column marginals are proportions that do not sum to 1 -
please respecify data")
    }
  }
  
  parties <- list()
  for (i in 1:(ncol(T))) {
    beta <- bayes.regress(T[,i] ~ G - 1, data = data, sample = sample,
                          weights = weights, truncate=truncate) 
    parties[[colnames(T)[i]]] <- beta
    colnames(parties[[colnames(T)[i]]]) <- colnames(G)
  }

 # ridx <- colnames(parties[[1]])
 # tmp.no <- matrix(NA, nrow(parties[[1]]), length(ridx))
 # colnames(tmp.no) <- ridx
 # for (i in 1:length(ridx)){
 #   tmp.sum <- rep(0, length(parties[[1]][,1]))
 #   for(j in 1:length(parties)){
 #     tmp.sum <- tmp.sum + parties[[j]][,i]
 #   }
 #   tmp.no[,ridx[i]] <- 1 - tmp.sum
 # }
 # NoVote <- tmp.no

  draws <- array(NA, dim=c(ncol(G), ncol(T), sample))
  for(i in 1:sample){
    for(j in 1:length(parties)){
      draws[,j,i] <- parties[[j]][i,]
    }
  #  draws[,(j+1),i] <- NoVote[i,]
  }
  rownames(draws) <- colnames(G)
  colnames(draws) <- colnames(T)

  out <- list(call=match.call(), draws=draws)
  class(out) <- "eiRegBayes"
  out
}
