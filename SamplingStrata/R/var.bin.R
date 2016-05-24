# ----------------------------------------------------
# Function to transform cotinuous stratification
# variables into ordered factors by using the
# k-means clustered method
# Author: Giulio Barcaroli
# Date: 4 January 2012
# ----------------------------------------------------
var.bin <- function(x, bins = 3, iter.max = 100) {
    km <- kmeans(x, bins, iter.max = 100)
    if (class(x) != "numeric" & class(x) != "integer") 
        stop("Kmeans clustering not applicable")
    km$cluster2 <- as.factor(km$cluster)
    levels(km$cluster2) <- levels(km$cluster2)[rank(km$centers)]
    var.binned <- as.numeric(as.character((km$cluster2)))
    return(var.binned)
}
