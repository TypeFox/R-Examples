# Popularity algorithm. (Recommandations ONLY)

mostpopular <- function(data) {
    
    x <- data@data
    
    colnames(x) <- NULL
    rownames(x) <- NULL
    
    items_nr_ratings <- colSums(x)
    names(items_nr_ratings) <- 1:ncol(x)
    items_nr_ratings <- sort(items_nr_ratings, decreasing = T)
    
    ordered_index <- as.numeric(names(items_nr_ratings))
    
    
    new("PPLclass", alg = "Popular", data = data, indices = ordered_index)
}

rrecsysRegistry$set_entry(alg = "Popular", 
                          fun = mostpopular, 
                          description = "Most popular algorithm.", 
                          reference = NA,
                          parameters = NA) 
