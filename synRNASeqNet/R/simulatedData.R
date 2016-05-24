simulatedData <-
function(p = 50, n = 100, mu = 100, sigma = 0.25, ppower = 1, noise = F, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  
  g <- barabasi.game(p, power = ppower, directed = F)
  #g <- erdos.renyi.game(n = p, p.or.m = 0.05, directed = FALSE)
  #plot(g)
  
  A <- as.matrix(get.adjacency(g, type = "both", names = T))
  diag(A) <- 1
  rownames(A) <- colnames(A) <- paste("gene", 1:p, sep = "_")
  
  nDistToCreate <- sum(A[upper.tri(A, diag = T)])
  muN <- rpois(nDistToCreate, mu)
  if(sigma == 0){
    poissGen <- t(sapply(muN, function(x) rpois(n, x)))
    #poissGen <- t(replicate(nDistToCreate, rpois(n, mu)))
  } else {
    rlognormpois <- function(n, lambda, sigma){
      eps <- rnorm(n)
      x <- rpois(n, exp(log(lambda) + sigma*eps))
      
      return(x)
    }
    
    nLognormpois <- ceiling(nDistToCreate * sample(seq(.10, .50, by = .1), size = 1))
    idx <- sample(1:nDistToCreate, size = nLognormpois)
    poissGen <- matrix(0, nrow = nDistToCreate, ncol = n)
    poissGen[-idx, ] <- t(sapply(muN[-idx], function(x) rpois(n, x)))
    poissGen[idx, ] <- t(sapply(muN[idx], function(x) rlognormpois(n, lambda = x, sigma = sigma)))
    #poissGen[-idx, ] <- t(replicate(nDistToCreate - nLognormpois, rpois(n, mu)))
    #poissGen[idx, ] <- t(replicate(nLognormpois, rlognormpois(n, lambda = mu, sigma = sigma)))
  }
  
  rowNames <- NULL
  for(i in 1:nrow(A)){
    tmp <- A[i, i:ncol(A), drop = F]
    rowNames <- c(rowNames, paste("G", i, "-", sub("gene_", "G", colnames(tmp[, which(tmp == 1), drop = F])), sep = ""))
  }
  rownames(poissGen) <- rowNames
  
  X <- matrix(0, p, n)
  rownames(X) <- rownames(A)
  colnames(X) <- paste("sample", 1:n, sep = "_")
  
  for(i in 1:p){
    tmp <- strsplit(rownames(poissGen), "-")
    tmp <- sapply(tmp, function(x) grep(paste("G", i, "$", sep = ""), x))
    tmp <- which(sapply(tmp, length) != 0)
    X[i, ] <- colSums(poissGen[tmp, , drop = F])
    if(noise == T){
      rrange <- max(X[i, ]) - min(X[i, ])
      X[i, ] <- X[i, ] + round(0.1 * (-1)^sample(0:1, size = n, replace = T) * sample(1:rrange, size = n, replace = T))
      X[i, X[i, ] <= 0] <- 0
    }
  }
  
  ans <- list(graph = g, adjMat = A, counts = X)
  return(ans)
}
