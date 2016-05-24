# Univariate:
EstimateIsingUni <- function(data, responses, beta = 1, ...){
  if (missing(responses)){
    responses <- sort(unique(c(data)))
  }
  
  if (length(responses) != 2){
    stop("Binary data required")
  }
  
  # Rescale data to binary:
  binarize <- function(x, responses){
    x2 <- x
    x2[x==responses[1]] <- 0
    x2[x==responses[2]] <- 1
    x2
  }
  data <- binarize(data, responses)
  
  # Number of variables:
  n <- ncol(data)
  
  # GLM for every node:
  Res <- lapply(seq_len(n), function(i) glm(data[,i] ~ data[,-i], family = binomial, ...))
  
  # Coefficients:
  Coefs <- lapply(Res, coef)
  
  # Thresholds:
  Thresholds <- sapply(Coefs, '[[', 1)
  # Network:
  Net <- matrix(0, n, n)
  for (i in seq_len(n)){
    Net[i,-i] <- Coefs[[i]][-1]
  }
  # Average:
  Net <- (Net + t(Net)) / 2
  
  # Rescale:
  Trans <- LinTransform(Net, Thresholds, c(0,1), responses)
  
  return(list(
    graph = Trans$graph,
    thresholds = Trans$thresholds,
    results = Res))
}




# Bivariate DOESNT WORK WHEN 11 COUNT IS LOW!:
EstimateIsingBi <- function(data, responses, beta = 1, ...){
  if (missing(responses)){
    responses <- sort(unique(c(data)))
  }
  
  if (length(responses) != 2){
    stop("Binary data required")
  }
  
  # Rescale data to binary:
  binarize <- function(x, responses){
    x2 <- x
    x2[x==responses[1]] <- 0
    x2[x==responses[2]] <- 1
    x2
  }
  data <- binarize(data, responses)
  
  # Number of variables:
  n <- ncol(data)
  
  # GLM for every pair of nodes:
  Res <- matrix(list(), n, n)
  Net <- matrix(0,n,n)
  ThresholdMat <- matrix(0,n-1,n)
  update <- rep(1, n)
  for (i in seq_len(n)){
    for (j in seq_len(i-1)){
      # Create factor:
      fac <- factor(data[,i] + 2*data[,j])
      Res[[i,j]] <- multinom(fac ~ data[,-c(i,j),drop=FALSE], ...)
      coef <- coef(Res[[i,j]])
      
      if (!is.matrix(coef)){
        coef <- as.matrix(t(coef))
        rownames(coef) <- unique(fac[fac!="0"])
      }
      
      # Add threshold i:
      if (1 %in% fac){
        ThresholdMat[update[i],i] <- coef[rownames(coef)=="1",1]
        update[i] <- update[i] + 1
      }
      
      # update threshold j:
      if (2 %in% fac){
        ThresholdMat[update[j],j] <- coef[rownames(coef)=="2",1]
        update[j] <- update[j]  + 1
      }
      
      # Update network:
      if (3 %in% fac){
        Net[i,j] <- Net[j,i] <- coef[rownames(coef)=="3",1]
      }
      
    }
  }
  
  # Threshold means:
  Thresholds <- colMeans(ThresholdMat)
  
  # Remove thresholds from net:
  Net <- Net - outer(Thresholds,Thresholds,'+')
  diag(Net) <- 0
  
  # Rescale:
  Trans <- LinTransform(Net, Thresholds, c(0,1), responses)
  
  return(list(
    graph = Trans$graph,
    thresholds = Trans$thresholds,
    results = Res))
}
# Contrasts matrix:
# Contr <- matrix(c(
#   -1, -1, 1,
#   1, -1, -1,
#   -1, 1, -1,
#   1, 1, 1  
# ),4,3, byrow = TRUE)
# rownames(Contr) <- 0:3
# colnames(Contr) <- 1:3
  
  
  
# As a loglinear model:
EstimateIsingLL <- function(data, responses, beta = 1, ...){
  if (missing(responses)){
    responses <- sort(unique(c(data)))
  }
  
  if (length(responses) != 2){
    stop("Binary data required")
  }
  
  # Rescale data to binary:
  binarize <- function(x, responses){
    x2 <- x
    x2[x==responses[1]] <- 0
    x2[x==responses[2]] <- 1
    x2
  }
  data <- binarize(data, responses)
  
  # Number of variables:
  n <- ncol(data)
  
  # GLM for every pair of nodes:
  data <- as.data.frame(data)
  for (i in seq_len(ncol(data))){
    data[[i]] <- factor(data[[i]], levels = c(1,0))
  }
  
  # Frequency table:
  Tab <- table(data)
  
  # margins:
  Margins <- alply(t(combn(seq_len(n),2)),1) 
  
  # Estimate
  Res <- loglin(Tab, Margins, param = TRUE, ...)
  
  # Parameters:
  Params <- Res$param
  Net <- matrix(0,n,n)
  Thresholds <- numeric(n)
  
  Names <- names(data)
  
  # Thresholds:
  for (i in seq_len(n)){
    Thresholds[[i]] <- Params[[Names[[i]]]][[1]]
  }
  
  # Network:
  for (j in seq_len(n)){
    for (i in seq_len(j-1)){
      Net[i,j]<- Net[j,i] <- Params[[paste0(Names[[i]],".",Names[[j]])]][1,1]
    }
  }
  
  # Rescale:
  Trans <- LinTransform(Net, Thresholds, c(-1,1), responses)
  
  return(list(
    graph = Trans$graph,
    thresholds = Trans$thresholds,
    results = Res))
}


