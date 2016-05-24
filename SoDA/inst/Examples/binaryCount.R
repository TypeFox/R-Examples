binaryCount <- function(nodes, leafValues) {
    nL <- length(leafValues)
    nN <- nrow(nodes)
    left <- nodes[,1]; right <- nodes[, 2]
    
    left <- ifelse(left<0, -left, left + nL)
    right <- ifelse(right<0, -right , right + nL)
    
    count <- c(leafValues, rep(NA, nN))
    
    while(any(is.na(count)))
        count <- c(leafValues, count[left] + count[right])

    count[-seq(length=nL)]
}
