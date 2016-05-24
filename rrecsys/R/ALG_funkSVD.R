# k is the number of features.
# gamma is the regularization term. 
# lambda is the learning rate.

# Reference: Y. Koren, R. Bell, and C. Volinsky. Matrix Factorization Techniques for Recommender Systems.
#            S. Funk. Netflix Update: Try this at Home.

FunkSVD <- function(data, k = 10, gamma = 0.015, lambda = 0.001) {
    
    x <- data@data
    
    colnames(x) <- NULL
    rownames(x) <- NULL
    
    row_x <- nrow(x)
    col_x <- ncol(x)
    if (ncol(x) < k | nrow(x) < k) 
        stop("Invalid number of features! \nLess features than the actual number of items or users! Please correct k!")
    # determine the means over row and columns and the total mean
    #biases <- calcBias(x)
    # initilize the user and item features
    U <- matrix(0.1, nrow = row_x, ncol = k)
    V <- matrix(0.1, nrow = col_x, ncol = k)
    
    #list of indices pointing to ratings on each item 
    itemIDX <- lapply(1:row_x, function(temp) which(x[temp, ] != 0))
    
    #list of indices pointing to ratings on each user 
    userIDX <- lapply(1:col_x, function(temp) which(x[, temp] != 0))
    
    # the training feature loop
    for (f in 1:k) {
        resetrrecsysenv()
        p <- matrix(100, nrow = row_x, ncol = col_x)
        # convergence check
        ptm <- Sys.time()
        while (!isConverged(x, p)) {

          p <- U %*% t(V) 
          error <- x - p
          
          # update user features
          
          temp_U <- U
          
          for (j in 1:col_x) {
            
            delta_Uik <- lambda * (error[userIDX[[j]], j] * V[j, f] - gamma * U[userIDX[[j]], f])
            
            U[userIDX[[j]], f] <- U[userIDX[[j]], f] + delta_Uik
            
          }
          
          # update item features
          for (i in 1:row_x) {
            
            delta_Vjk <- lambda * (error[i, itemIDX[[i]]] * temp_U[i, f] - gamma * V[itemIDX[[i]], f])
            
            V[itemIDX[[i]], f] <- V[itemIDX[[i]], f] + delta_Vjk
            
          }

        }  #end convergence loop
        
        cat("Feature:", f, "/", k, " trained. Time:", as.numeric(Sys.time() - ptm, units = "secs"), "seconds. \n")
        
    }  #end feature loop
    
    p_FunkSVD <- list(k = k, gamma = gamma, lambda = lambda)
    
    new("SVDclass", alg = "FunkSVD", data = data, factors = list(U = U, V = V), parameters = p_FunkSVD)
    
}


p_FunkSVD <- list(k = 10, gamma = 0.01, lambda = 0.001)
rrecsysRegistry$set_entry(alg = "FunkSVD", 
                          fun = FunkSVD, 
                          description = "Funk SVD", 
                          reference =  "Y. Koren, R. Bell, and C. Volinsky. Matrix Factorization Techniques for Recommender Systems. \nS. Funk. Netflix Update: Try this at Home.",
                          parameters = p_FunkSVD) 

#           INTERACTIVE IMPLEMANTATION OF THE UPDATES ON THE FEATURES
#             for (i in 1:nrow(x)){
#               
#               temp_Ui <- U[i,]
#               
#               for(j in itemIDX[[i]]){
#                 eij <- x[i,j] - sum(U[i,] * V[j,])
#                 delta_Uik <- lambda * (eij * V[j, f] - gamma * U[i,f])
#                 U[i, f] <- U[i, f] + delta_Uik
#                 delta_Vjk <- lambda * (eij * temp_Ui[f] - gamma * V[j , f])
#                 V[j , f] <- V[j , f] + delta_Vjk
#                 
#               }
#             }
#           p <- U %*% t(V)
