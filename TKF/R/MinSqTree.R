validMatrix <- function(d){
  stopifnot(nrow(d) == ncol(d))
  if(is.null(colnames(d))){
    colnames(d) <- 1:nrow(d)
  }
  if(is.null(rownames(d))){
    rownames(d) <- 1:nrow(d)
  }
  stopifnot(identical(colnames(d), rownames(d)))
  return(d)
}

## function performs weighted least-squares phylogeny inference by nni
optim.phylo.wls <- function(Dist, Var=NULL, stree=NULL, set.neg.to.zero=TRUE,
                           fixed=FALSE, tol=1e-10, collapse=TRUE){
  ## ordinary unweighted least square tree
  if(is.null(Var)){
    ans <- optim.phylo.ls(Dist, stree=stree, set.neg.to.zero=set.neg.to.zero,
                          fixed=fixed, tol=tol, collapse=collapse)
    return(ans)
  }

  ## weighted least square tree
  # change Dist and Var to a matrix (if actually an object of class "dist")
  if(is(Dist, "dist")){
    Dist <- as.matrix(Dist)
  }
  if(is(Var, "dist")){
    Var <- as.matrix(Var)
  }

  # Dist and Var must have same dimetions.
  stopifnot(all(dim(Dist) == dim(Var)))
  Dist <- validMatrix(Dist)
  Var <- validMatrix(Var)


  # compute the number of species
  n<-nrow(Dist)

  if(!is(stree, "phylo")){
    message("random starting tree if stree is not phylo class!")
    stree <- rtree(n=n, tip.label=rownames(Dist), br=NULL, rooted=FALSE)
  }

  if(!is.binary.tree(stree)){
    stree <- multi2di(stree)
  }
  if(is.rooted(stree)){
    stree <- unroot(stree)
  }
  
  # get ls branch lengths for stree
  best.tree <- ls.tree(stree, Dist, Var)
  Q <- attr(best.tree, "Q-score")
  bestQ <- 0 # to start the loop

  # for search
  Nnni <- 0

  # loop while Q is not improved by nni
  while(bestQ - Q < tol && fixed==FALSE){

    nni.trees <- nni(best.tree)
    nniQ <- vector()

    bestQ <- Inf
    for(i in 1:length(nni.trees)){

      # compute least squares branch lengths and Q
      nni.trees[[i]] <- ls.tree(nni.trees[[i]], Dist, Var)

      # compute Q
      nniQ[i] <- attr(nni.trees[[i]], "Q-score")

      # is this the best one so far?
      if(nniQ[i] < bestQ){
        bestQ <- nniQ[i]
        ind <- i
      }
    }

    # set new best tree
    if(bestQ < Q){
      best.tree <- nni.trees[[ind]]
      Nnni <- Nnni + 1
      Q <- attr(best.tree, "Q-score")
      print(paste(Nnni,"set(s) of nearest neighbor interchanges. best Q so far =",round(Q,10),collapse=""))
    } else bestQ <- Inf
  }
  message(paste("best Q score of",round(Q,10),"found after",Nnni,"nearest neighber interchange(s).",collapse=""))
  if(set.neg.to.zero){
    best.tree$edge.length[best.tree$edge.length < 0] <- 0
  }

  attr(best.tree, "Q-score") <- Q
  if(collapse){
    best.tree <- di2multi(best.tree)
  }
  return(best.tree)

}

# function computes the ls branch lengths and Q score for a tree
# written by Liam J. Revell 2011

ls.tree <- function(tree, Dist, Var){
  Dist <- validMatrix(Dist)
  Var <- validMatrix(Var)
  stopifnot(identical(dimnames(Dist), dimnames(Var)))
  # We only need to check one rownames of Dist or Var.
  stopifnot(setequal(tree$tip.label, rownames(Dist)))
  
  # compute design matrix for tree i
  X <- phyloDesign(tree)
  
  # sort and columnarize Dist and Var
  Dist <- Dist[tree$tip.label, tree$tip.label]
  colDist <- Dist[lower.tri(Dist)]
  Var <- Var[tree$tip.label, tree$tip.label]
  colVar <- Var[lower.tri(Var)]

  # compute the least squares branches conditioned on tree i
  #v <- solve(t(X)%*%X) %*% t(X) %*% colDist
  toFit <- as.data.frame(cbind(X, colDist))
  fit <- lm(colDist~0+., data=toFit, weights=1/colVar)
  v <- coef(fit)
  v <- matrix(v, dimnames=list(colnames(X)), ncol=1) 
  tree$edge.length <- v

  # compute the distances for this tree
  d <- X %*% v

  # compute Q
  Q <- sum((fit$residuals)^2 * fit$weights)

  # assign attribute to tree
  attr(tree, "Q-score") <- Q

  return(tree)
}

