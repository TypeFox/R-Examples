

#' Local Indicator of Phylogenetic Association
#'
#' This function computes Local Indicator of Phylogenetic Association (local Moran's I)
#' for each tip of a tree.
#' Tests are based on permutations.
#' 
#' @param p4d a \code{phylo4d} object.
#' @param trait the traits in the \code{phylo4d} object for which to compute LIPA.
#' Can be a character vector giving the name of the traits
#' or numbers giving the column in the table of the data slot of the \code{phylo4d} object.
#' @param reps a numeric value. Number of repetitions for the estimation of p.values with randomization.
#' @param alternative a character string specifying the alternative hypothesis for the tests.
#' Must be one of \code{greater} (default), \code{two-sided} or \code{less}.
#' @param prox.phylo a matrix of phylogenetic proximities or
#' a character string specifying a method to compute it.
#' See Details.
#' @param as.p4d logical. Should the results returned as a \code{phylo4d} object?
#' 
#' @details If "\code{prox.phylo}" is a character string,
#' the phylogenetic proximity matrix is computed internally
#' using the function \code{\link[adephylo]{proxTips}} from the package \pkg{adephylo}.
#' Different methods are available:
#' "\code{patristic}","\code{nNodes}","\code{Abouheif}" and "\code{sumDD}".
#' See \code{\link[adephylo]{proxTips}} for details about the methods.
#' 
#' @return If \code{as.p4d} is \code{FALSE} (default), the function returns a list:
#' \describe{
#'   \item{lipa}{A matrix of LIPA indices computed for each tip of the tree and each trait.}
#'   \item{p.value}{A matrix of p-values (tests of LIPA indices)}
#'   \item{reps}{Number of permutations for the tests}
#'   \item{alternative}{Alternative hypothesis for the tests}
#' }
#' If \code{as.p4d} is \code{TRUE}, the function returns a \code{phylo4d} object
#' with LIPA values as tips associated data.
#'
#' @references Anselin L. (1995) Local Indicators of Spatial Association-LISA. Geographical Analysis 27, 93-115.

#'
#' @export
lipaMoran <- function(p4d, trait = names(tdata(p4d)), reps=999,
                 alternative = "greater", prox.phylo = "patristic", as.p4d = FALSE){
  
  p4 <- extractTree(p4d)
  phy <- as(p4, "phylo")
  new.order <- phy$edge[, 2][!phy$edge[, 2] %in% phy$edge[, 1]]
  tips <- phy$tip.label[new.order]
  n.tips <- length(tips)
  X <- tdata(p4d, type = "tip")
  X <- X[tips, trait]
  X <- as.data.frame(X)
  colnames(X) <- trait
  n.traits <- ncol(X)
  if(is.numeric(trait)){
    trait <- names(tdata(p4d))[trait]
  }
  if(is.vector(prox.phylo) & is.character(prox.phylo)){
    prox.phylo <- match.arg(prox.phylo, c("patristic", "nNodes", "Abouheif", "sumDD"))
    W <- proxTips(phy, method = prox.phylo)[tips, tips]
  } else {
    if(is.matrix(prox.phylo)){
      W <- prox.phylo[tips, tips]
    } else {
      stop("prox.phylo is not valid")
    }
  }
  X <- scale(X, scale = FALSE)
  X.lag <- W %*% X
  S2 <- diag(var(X))
  rt <- t(X) / S2
  lipa <- X.lag * t(rt)
  
  if(reps > 0){
    lipa.perm <- array(NA, dim=c(dim(X), reps))
    for(i in 1:reps){
      Xs <- apply(X, 2, sample)
      Xs.lag <- W %*% Xs
      rts <- t(Xs) / S2
      lipa.perm[ , , i] <- Xs.lag * t(rts)
    }
    
    alternative <- match.arg(alternative, c("two-sided", "greater", "less"))
    p.value <- matrix(ncol = n.traits, nrow = n.tips, dimnames = list(tips, trait))

    for(i in 1:n.tips){
      for(j in 1:n.traits){
        
        if(alternative == "two-sided"){
          lipa.perm.mean <- mean(lipa.perm[i, j, ])
          lipa.perm0 <- abs(lipa.perm[i, j, ] - lipa.perm.mean)
          lipa.obs0 <- abs(lipa[i, j] - lipa.perm.mean)
          p.value[i, j] <- (sum(lipa.perm0 >= lipa.obs0) + 1)/(reps + 1)
        }
        
        if(alternative == "greater"){
          p.value[i, j] <- (sum(lipa.perm[i, j, ] >= lipa[i, j]) + 1)/(reps + 1)
        }
        
        if(alternative == "less"){
          p.value[i, j] <- (sum(lipa.perm[i, j, ] <= lipa[i, j]) + 1)/(reps + 1)
        }
      }
    }
  }

  if(as.p4d){
    res <- phylo4d(phy, lipa)
  } else {
    if(reps > 0){
      res <- list(lipa = lipa, p.value = p.value, reps = reps, alternative = alternative)
    } else {
      res <- lipa
    }
  }
  
  return(res)
}

