`dist.eval` <-
function(x, dist){
    x <- as.matrix(x)
    tots <- rowSums(x)
    if(any(tots==0)) {
        cat("Warning: the community matrix contains some sites with only zero abundances\n")
        cat("You may want to use functions removezerospecies or dist.zeroes from Biodiversity.R\n")
    }else{
        tests <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", 
            "chao", "cao", "mahalanobis")
        op <- options()
        options(warn=-1)
        if (any(dist==tests)) {
            dist1 <- vegdist(x, method=dist)
        }else{
            if((dist %in% c("w", "-1", "c", "wb", "r", "I", "e", "t", "me", "j", "sor", "m", "-2", "co", "cc", "g", "-3", "l", "19", "hk", "rlb", "sim", "gl", "z")) == F) {
                stop("Provide acceptable method for betadiver")
            }
            dist1 <- betadiver(x, method=dist)
        }
        dist2 <- no.shared(x)
        list1 <- (dist2==0)
        list2 <- (dist2==1)
        max <- max(dist1[list1])
        min <- min(dist1[list2])
        options(op)
        if(min<max) {
            cat("Warning: min distance for sites with no shared species(",min,") < max dist for other sites(", max, ")\n")
            cat("Choose other distance measure or stepacross\n")
        }
        return(distconnected(dist1))
    }
}

