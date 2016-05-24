recommend <- function(model, topN = 3) {
  
    x <- model@data@data
        
    if (topN >= ncol(x)) 
        stop("topN value is larger than the number of items that can be recommended.")
    if (topN < 1) 
        stop("Not valid value for topN.")
    
    rated_items <- which(x != 0)
    
    item_names <- colnames(x)
    
    if (model@alg == "Popular") {

        rec_indices <- lapply(1:nrow(x), function(m) model@indices[1:topN])
        
        if (is.null(item_names)) {
            recnames <- rec_indices
        } else {
            recnames <- lapply(1:nrow(x), function(i) item_names[rec_indices[[i]]])
        }
        names(recnames) <- rownames(x)
        
        
        
        return(new("recResultsClass", indices = rec_indices, recommended = recnames))
    }
    
    p <- predict(model)
    
    colnames(p) <- NULL
    
    p[rated_items] <- NA
    
    rec_indices <- lapply(1:nrow(x), function(i) order(p[i, ], na.last = NA, decreasing = T)[1:topN])
    rec_indices <- lapply(1:nrow(x), function(i) rec_indices[[i]][!is.na(rec_indices[[i]])])
    
    if (is.null(item_names)) {
        recnames <- rec_indices
    } else {
        recnames <- lapply(1:nrow(x), function(i) item_names[rec_indices[[i]]])
    }
    names(recnames) <- rownames(x)
    
    new("recResultsClass", indices = rec_indices, recommended = recnames)
    
} 


#rated_items <- lapply(1:nrow(x), function(i) which(x[i, ] != 0))

#rec_indices <- lapply(rated_items, function(m) (model@indices[!model@indices %in% m])[1:topN])
