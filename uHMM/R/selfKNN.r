# KNN dans un meme jeu de donnees
#' @title Self KNN 
#' @description This function performs the k-Nearest Neighbour algorithm without class estimation, but only computation of distances and neighbours.
#' @param train numeric matrix or data frame. 
#' @param K number of neighbours considered.
#' @return The function returns a list with the following components:
#' \item{D}{matrix of squared root of the distances between observations and their nearest neighbours.}
#' \item{idx}{Index of K nearest neighbours of each observation.}
#' @export
#' @importFrom stats quantile
#' @examples 
#' x<-matrix(runif(10),ncol=2)
#' plot(x,pch=c("1","2","3","4","5"))
#' selfKNN(x,K=4)

selfKNN <- function(train, K=1){
  
  idx <- matrix(0, nrow(train), K)
  D <- idx
    
    dst=as.matrix(dist(train,upper=TRUE))
    diag(dst)<-Inf
    n<-nrow(train)
    
    for(k in 1:nrow(train)){
      obs<-dst[,k]
      
      if (K==1){
        idx[k] = which.min(obs)
        D[k] = obs[idx[k]]
      }else{
        #tri partiel moins long que tri complet
        quant <- quantile(obs, probs=(K+1)/n) #+1 pour etre sur d'en avoir K
        quant.idx <- which(obs<=quant)
        idx[k,] <- quant.idx[order(obs[quant.idx], decreasing=FALSE)[1:K]]
        D[k,] <- obs[idx[k,]]
      }
    }

  return(list(D=D, idx=idx))
}