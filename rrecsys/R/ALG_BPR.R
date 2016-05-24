# lambda: the learning rate. 
# regU: regularization parameter for user factors. 
# regI: regularization parameter for positive item factors. 
# regJ: regularization parameter for negative item factors. 
# updateJ: update factors for negative sampled items during learning(default value TRUE).

# Reference: S. Rendle, C. Freudenthaler, Z. Gantner, and L. Schmidt-Thieme. BPR: Bayesian Personalized Ranking from Implicit Feedback.

BPR <- function(data, k = 10, lambda = 0.05, regU = 0.0025, regI = 0.0025, regJ = 0.0025, updateJ = TRUE) {
    
    x <- data@data
    
    row_x <- nrow(x)
    col_x <- ncol(x)
    
    colnames(x) <- NULL
    rownames(x) <- NULL
    
    if (ncol(x) < k | nrow(x) < k) 
        stop("Invalid number of features! \nLess features than the actual number of items or users! Please correct k!")
    
    # initilize the user and item features
    U <- matrix(0.1, nrow = row_x, ncol = k)
    V <- matrix(0.1, nrow = col_x, ncol = k)
    
    #list of indices pointing to ratings on each user 
    userIDX <- lapply(1:row_x, function(i) which(x[i, ] >= data@minimum))
    userIDX <- lapply(userIDX, unname)
    #list of indices pointing to unrated items on each user 
    userNOIDX <- lapply(1:row_x, function(i) which(x[i, ] < data@minimum))
    userNOIDX <- lapply(userNOIDX, unname)
    
    p <- U %*% t(V)
    resetrrecsysenv()
    
    while (!isConverged(x, p)) {
        
        for (s in 1:100 * row_x) {
            # extract a random user one random rated item and one random unrated item for that user.
            while (TRUE) {
                u <- sample(1:row_x, 1)
                if (length(userIDX[[u]]) == 0) 
                  next
                i <- userIDX[[u]][sample(1:length(userIDX[[u]]), 1)]
                j <- userNOIDX[[u]][sample(1:length(userNOIDX[[u]]), 1)]
                break
            }
            
            # predict xui and xuj
            xui <- sum(U[u, ] * V[i, ])
            xuj <- sum(U[u, ] * V[j, ])
            
            xuij <- xui - xuj
            
            sigma <- 1/(1 + exp(-xuij))
            
            loss <- -log(1/(1 + exp(xuij)))
            
            U[u, ] <- U[u, ] + lambda * (sigma * (V[i, ] - V[j, ]) - regU * U[u, ])
            
            V[i, ] <- V[i, ] + lambda * (sigma * U[u, ] - regI * V[i, ])
            
            if (updateJ) {
                V[j, ] <- V[j, ] + lambda * (sigma * (-U[u, ]) - regJ * V[j, ])
            }
            
            
        }
        
        p <- U %*% t(V)
        
    }  #convergence
    
    p_BPR <- list(k = k, lambda = lambda, regU = regU, regI = regI, regJ = regJ, updateJ = updateJ)
    
    new("BPRclass", alg = "BPR", data = data, factors = list(U = U, V = V), parameters = p_BPR)
}


p_BPR <- list(k = 10, lambda = 0.05, regU = 0.0025, regI = 0.0025, regJ = 0.0025, updateJ = TRUE)
rrecsysRegistry$set_entry(alg = "BPR", 
                          fun = BPR, 
                          description = "Bayesian Personalized Ranking.", 
                          reference = "S. Rendle, C. Freudenthaler, Z. Gantner, and L. Schmidt-Thieme. BPR: Bayesian Personalized Ranking from Implicit Feedback.",
                          parameters = p_BPR) 
