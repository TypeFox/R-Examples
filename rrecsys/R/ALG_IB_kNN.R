# neigh: neighbourhood size.

# Reference: B. Sarwar, G. Karypis, J. Konstan, and J. Riedl. Item-based collaborative filtering recommendation algorithms.

IB_kNN <- function(data, neigh = 10) {

    maxNumRatings <- max(rowRatings(data))
    x <- data@data
    
    colnames(x) <- NULL
    rownames(x) <- NULL
    
    ptm <- Sys.time()
    
    if (neigh > ncol(x) - 1) 
        stop("Invalid value for k!!! Please change the k attribute.\nNeighborhood value is larger than the entire number of items in the dataset.")
    
    if (neigh < 1) 
        stop("Invalid value for neigh!!!")
    
    rated_index <- lapply(1:nrow(x), function(q) which(x[q, ] != 0))

    n_x <- x 

    temp <- rowRatings(data)
    
    means <-  rowSums(x)/temp    #rowMeans(x) #apply(x, 1 ,function(m) sum(m)/maxNumRatings) #
    for(i in 1:nrow(x)) n_x[i, rated_index[[i]]] <- n_x[i,rated_index[[i]]] - means[i] 
    sim <- simil(n_x, method = 'cosine',diag = TRUE, by_rows = FALSE) 
    sim <- as.matrix(sim) 
    colnames(sim) <- NULL 
    rownames(sim)<- NULL
    sim_index_kNN <- t(apply(sim, 1, function(q) 
      order(q, decreasing = TRUE, na.last = TRUE)))[, 1:neigh]
    
    # every item is similar 100% to itself. 
    diag(sim) <- 1
    
    cat("Neighborhood calculated in: ", as.numeric(Sys.time() - ptm, units = "secs"), "seconds.\n")
    
    if (neigh == 1) 
        sim_index_kNN <- as.matrix(sim_index_kNN)
    
    new("IBclass", alg = "IBKNN", data = data, sim = sim, sim_index_kNN = sim_index_kNN, neigh = neigh)
    
}


rrecsysRegistry$set_entry(alg = "IBKNN", 
                          fun = IB_kNN, 
                          description = "Item based k-NN", 
                          reference = "B. Sarwar, G. Karypis, J. Konstan, and J. Riedl. Item-based collaborative filtering recommendation algorithms.",
                          parameters = list(neigh = 10)) 




#     # adjusted cosine similarity
#     sim <- matrix(NA, nrow = ncol(x), ncol = ncol(x))
#     rowM <- apply(x, 1, function(m) sum(m)/maxNumRatings)
#     
#     n_x <- x
#     for (i in 1:nrow(x)) n_x[i, rated_index[[i]]] <- n_x[i, rated_index[[i]]] - rowM[i]
#     
#     # simility matrix is symmetric so only the lower triagnle needs to be computed than it can be projected to the other half
#     for (i in 2:ncol(x)) {
#         for (j in 1:(i - 1)) {
#           
#           corated <- which((n_x[, i] != 0) & (n_x[, j] != 0))
#           
#           if(length(corated) >= min_corated){
#             numerator <- crossprod(n_x[corated, i], n_x[corated, j])
#             denominator <- sqrt(crossprod(n_x[corated, i]) * crossprod(n_x[corated, j]))
#             if (denominator != 0) {
#               sim[i, j] <- numerator/denominator
#             }
#             sim[j, i] <- sim[i, j]
#           }
#         }
#     }
#   
#     sim_index_kNN <- t(apply(sim, 1, function(q) order(q, decreasing = TRUE, na.last = TRUE)))[, 1:neigh]
#     
