
computeDelta <- function(X, ExpInd, covariance=TRUE){
  ## difference between the covariance matrices of X under the different interventions

    n <- nrow(X)
    settings <- sort(unique(ExpInd))
    G <- length(settings)
    
    # initialize list containing Delta matrices
    Deltalist <- vector("list", G)
    
    for (k in 1:G){
        useo <- which(ExpInd==settings[k])
        useo2 <- which(ExpInd!=settings[k])
        if(covariance){
            Deltalist[[k]] <- cov( X[useo ,]) - cov( X[useo2,])
        }else{
            Deltalist[[k]] <- 1/(length(useo))*(t(X[useo,])%*%X[useo,]) - 1/(length(useo2))*(t(X[useo2,])%*%X[useo2,])
        }
    }
  
    list(Delta = Deltalist, K = G)
}

