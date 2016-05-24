#Algorithm that generates predictions based on the global average only.


globalAverage <- function(data) {
    
    average <- sum(data@data) / numRatings(data)
    
    average <- matrix(1, nrow = nrow(data), ncol = ncol(data)) * average
    
    new("algAverageClass", alg = "globalAverage", data = data, average = average)
    
}

#Algorithm that generates predictions based on the items average only.
itemAverage <- function(data) {
    
  
    average <- colSums(data@data) / colRatings(data)

    names(average) <- NULL
    
    if (any(is.nan(average))) {
        average[which(is.nan(average))] <- sum(data@data) / numRatings(data)
    }
    
    average <- matrix(rep(average, nrow(data)), nrow = nrow(data), byrow = TRUE)
    
    new("algAverageClass", alg = "itemlAverage", data = data, average = average)
    
}

#Algorithm that generates predictions based on the user average only.
userAverage <- function(data) {
    
    average <- rowSums(data@data) / rowRatings(data)
    
    names(average) <- NULL
    
    if (any(is.nan(average))) {
        average[which(is.nan(average))] <- sum(data@data) / numRatings(data)
    }
    
    average <- matrix(rep(average, ncol(data)), nrow = nrow(data), byrow = F)
    
    
    new("algAverageClass", alg = "userAverage", data = data, average = average)
    
}

rrecsysRegistry$set_entry(alg = "itemAverage", 
                          fun = itemAverage, 
                          description = "Item average",
                          reference = NA,
                          parameters = NA)

rrecsysRegistry$set_entry(alg = "userAverage", 
                          fun = userAverage, 
                          description = "User average",
                          reference = NA,
                          parameters = NA)

rrecsysRegistry$set_entry(alg = "globalAverage", 
                          fun = globalAverage, 
                          description = "Global average", 
                          reference = NA,
                          parameters = NA) 
