defineData <- function(data, binary = FALSE, minimum = 1, maximum = 5, halfStar = FALSE, goodRating = 0.5) {
    
    data[is.nan(data)] <- 0
    data[is.null(data)] <- 0
    
    if (!binary) {
        new("dataSet", data = data, binary = binary, minimum = minimum, maximum = maximum, halfStar = halfStar)
    } else {
        data[data < goodRating] <- 0
        data[data >= goodRating] <- 1
        new("dataSet", data = data, binary = binary, minimum = 1, maximum = 1, halfStar = FALSE)
    }
} 

#dataSet####
setMethod("show", signature(object = "dataSet"), function(object) {
  if (!object@binary) {
    cat("Dataset ")
  } else {
    cat("Binary dataset ")
  }
  cat("containing", nrow(object@data), "users and ", ncol(object@data), "items.")
})