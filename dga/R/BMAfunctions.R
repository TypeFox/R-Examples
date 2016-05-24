#     DGA does capture-recapture estimation using decomposable graphical models
#
#     Copyright (C) 2014, Human Rights Data Analysis Group (HRDAG)
#     https://hrdag.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#library(R.utils)
#library(network)
#functions to do the madigan + york bma routine


CompLogML <- function(D, delta){
  #function to calculate the log marginal likelihood of a graph component
  #operates on rows
  #need to test this against the matlab code still
  #D should be the "squished" marginal table
  #delta is the prior cell counts
  #if (dim(D)[2]%%2 != 0){ print('problem! This is not a marginal table!')}
  out <- apply(lgamma(D + delta), 1, sum) - apply(lgamma(D*0 + delta), 1, sum)
  return(out)
}

integer.base.b <- function(x, b=2){
    xi <- as.integer(x)
    if(any(is.na(xi) | ((x-xi)!=0)))
      print(list(ERROR="x not integer", x=x))
    N <- length(x)
    xMax <- max(x)
    ndigits <- (floor(logb(xMax, base=2))+1)
    Base.b <- array(NA, dim=c(N, ndigits))
    for(i in 1:ndigits){#i <- 1
      Base.b[, ndigits-i+1] <- (x %% b)
      x <- (x %/% b)
    }
    if(N ==1) Base.b[1, ] else Base.b
  }

MakeCompMatrix <- function(p, delta, Y, Nmissing){
  #pre-computes each component's  LogML
  compLMLs <- matrix(0,nrow = 2^p-1, ncol = length(Nmissing))
  bins <- integer.base.b(1:(2^p-1), 2)
  for(i in 1:(2^p-1)){
    inds <- which(bins[i,]==1)
    D <- c(apply(Y, inds, sum))
    Dmat <- t(matrix(D, ncol = length(Nmissing), nrow = length(D)))
    Dmat[,1] <- Dmat[,1] + Nmissing

    #compute alpha
    alpha <- rep(delta * 2^(p - sum(bins[i,])), ncol(Dmat))
    compLMLs[i,] <- CompLogML(Dmat, alpha)
  }
  return(compLMLs)
}

tmpfun <- function(x,p){
  tmp <- rep(0,p)
  tmp[x] <- 1
  return((tmp))
}


bma.cr <- function(Y, Nmissing, delta, graphs, logprior=NULL, log.prior.model.weights=NULL, normalize=TRUE){
  #function to madigan + york style bma
  #Y is the array of counts.. should be dimension 2x2x2 (if p = 3)
  #Nmissing is the set of possible uncounted cases
  #delta is the prior weight for each cell
  #graphs are all of the decomposable graphs for p lists. These are pre-computed for p = 3, 4, 5.
  #normalize indicates whether to normalize the count/model weights prior to output. Default is to normalize, this is mostly used for development.
  #log.prior.model.weights: a vector of prior model weights. Should be the same length as length(graphs)

  #get number of lists
  p <- length(dim(Y))
  
  #set the missing cell to be 0
  Y[1] <- 0 
  
  #model x estimate weights go in here
  modNweights <- matrix(nrow = length(graphs), ncol = length(Nmissing))
  
  
  #first pre-compute the matrix of component-wise LMLs
  compMat <- MakeCompMatrix(p, delta, Y, Nmissing) # all but the last graph (the one that doesn't really matter) match with matlab code for 3 lists

  j <- 1
  for(graph in graphs){#loop over all possible models
    #graph$C cliques of the graph
    #graph$S separators
    binC <-t(sapply(graph$C, tmpfun, p = p))
    decC <- apply(t(binC)*rev(2^(0:(p-1))), 2, sum)
    compMats <- compMat[decC,]
    if(!is.null(nrow(compMats))){
      cliqueML <- apply(compMats, 2, sum)
    }else{cliqueML <- compMats}

    if(!is.null(graph$S)){
      binS <-t(sapply(graph$S, tmpfun, p = p))
      decS <- apply(t(binS)*rev(2^(0:(p-1))), 2, sum)

      compMats <- compMat[decS,]
      if(!is.null(nrow(compMats))){
        sepML <- apply(compMats, 2, sum)
      }else{sepML <- compMats}
    }else{sepML <- 0; decS <- NULL}

    nsubgraphs <- length(decC) - length(decS)
    #nsubgraph.*(gammaln(sum(alpha))-gammaln(N+sum(alpha)));

    alpha <- rep(delta, length(Y))
    modNweights[j,] <- (cliqueML - sepML + nsubgraphs*(lgamma(sum(alpha)) - lgamma(Nmissing + sum(Y) + sum(alpha))))
    j <- j + 1
  }


  #add on multinomial coefficient
  multicoef <- lgamma(Nmissing + sum(Y) + 1) - sum(lgamma(Y[-1] + 1)) - lgamma(Nmissing + 1)
  #multicoef <- lgamma(Nmissing + sum(Y) ) - sum(lgamma(Y[-1] )) - lgamma(Nmissing )

  modNweights <- t(t(modNweights) + multicoef)

  #add on prior
  if(is.null(logprior)){
  logprior <- -log(sum(Y) + Nmissing)}

  modNweights <- t(t(modNweights) + logprior) 

  if(!is.null(log.prior.model.weights)){
    modNweights <- modNweights + log.prior.model.weights
  }

  if(normalize){
    modNweights <- modNweights - max(modNweights)
    weights <- exp(modNweights)
    weights <- weights/sum(weights)
  } else{ 
      weights <- modNweights
    }
  return(weights)
}

plotPosteriorN <- function(weights, N, main=NULL){
  #this function
  plot(N, apply(weights, 2, sum), type = 'l', col = 'black', lwd = 3, ylab = "Posterior Probability of N", xlab = "N", ylim=c(0, 1.25*max(apply(weights, 2, sum))))
  title(main)
  wts <- apply(weights, 1, sum)
  for(i in 1:nrow(weights)){
    lines(N, weights[i,], lwd = wts[i]*3, lty = 'dashed')
  }
  legend("topright", legend = c("Averaged Post. Prob.", "Post. Prob. By Model"), lty  = c(1, 2), cex = .75)
}

#plotTopModels <- function(graphs, weights, p, how.many = 3){
#  modelWeights <- apply(weights, 1, sum)
#  or <- rev(order(modelWeights))
#  par(mfrow=c(how.many,1))
#  for(i in 1:how.many){
#    graph <- graphs[[or[i]]]
#    Adj <- makeAdjMatrix(graph,p)
#    rownames(Adj) <- c("one", "two", "three")
#    net <- as.network(Adj)
#    plot(net, main = paste("posterior prob ", round(modelWeights[or[i]], 3)))
#    #still doesn't plot names, need to figure out how to get it to do this
#  }
#}

makeAdjMatrix <- function(graph, p){
  Adj <- matrix(0, nrow = p, ncol = p)
  diag(Adj) <- 1
  for(i in 1:length(graph$C)){
    if(length(graph$C[[i]])>1){
      combns <- combn(graph$C[[i]],2)
      Adj[combns[1], combns[2]] <- 1
    }
  }
  for(i in 1:length(graph$S)){
    if(length(graph$S[[i]])>1){

      combns <- combn(graph$S[[i]],2)
      Adj[combns[1], combns[2]] <- 1
    }
  }

  Adj <- Adj + t(Adj)
  Adj[Adj>1] <- 1
  return(Adj)
}

