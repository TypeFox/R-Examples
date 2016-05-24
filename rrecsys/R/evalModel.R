setMethod("evalModel", signature = c(data = "dataSet"), function(data, folds) {
    
    if (missing(folds)) 
        folds <- 5
    
    if (folds < 2) 
        stop("k-fold cross validation requires at least two folds to continue!")
    
    x <- data@data
    
    fold_indices <- vector("list", folds)
    
    fold_indices_x_user <- vector("list", folds)
    
    fold_indices_x_user <- lapply(fold_indices_x_user, function(q) q <- vector("list", nrow(x)))
    
    item_names <- colnames(x)
    
    colnames(x) <- NULL
    
    rated_index_by_row <- lapply(1:nrow(x), function(temp) which(x[temp, ] != 0))
    
    if (folds <= 1) 
        stop("Invalid number of folds!!! k - fold attirubute is lower than 0.")
    
    if (folds > length(which(x != 0))) 
        stop("Invalid number of folds!!! k is bigger than the number of ratings. ")
    
    fold_vec <- c(1:folds)
    
    for (i in 1:length(rated_index_by_row)) {
        
        while (length(rated_index_by_row[[i]]) > 0) {
            
            # select a random fold
            a_fold <- sample(1:length(fold_vec), 1)
            j <- fold_vec[a_fold]
            fold_vec <- fold_vec[-a_fold]
            
            # choosing the item to put in the fold and extract its index
            an_item <- sample(1:length(rated_index_by_row[[i]]), 1)
            temp <- rated_index_by_row[[i]][an_item]
            rated_index_by_row[[i]] <- rated_index_by_row[[i]][-an_item]
            
            fold_indices_x_user[[j]][[i]] <- c(temp, fold_indices_x_user[[j]][[i]])
            
            # storing the index in the folds
            temp <- (temp - 1) * nrow(x) + i
            
            fold_indices[[j]] <- c(temp, fold_indices[[j]])
            
            if (length(fold_vec) == 0) {
                fold_vec <- c(1:folds)
            }
        }
    }
    
    colnames(x) <- item_names
    
    new("evalModel", data = data, folds = folds, fold_indices = fold_indices, fold_indices_x_user = fold_indices_x_user)
}) 
