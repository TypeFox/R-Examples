MutualInformation <-
function(Q, px = NA) {
    ## This function computes the mutual information for a given
    ## channel, optionally using a specified source distribution.
    
    m <- nrow(Q$x);
    n <- nrow(Q$y);
    R <- 0.0;
    
    if(is.na(px[1])) {
        px <- Q$px
    } else {
        if(length(px) != m) {
            stop("Alphabet for channel and source must match.")
        }
    }

    cpd <- vector(mode = "list", length = m)
    for(i in 1:m) {
        cpd[[i]] <- ConditionalDistribution(Q, i)
    }
    
    qy <- vector(mode = "numeric", length = n)
    for(j in 1:n) {
        qy[j] <- 0
        for(i in 1:m) {
            qy[j] <- qy[j] + cpd[[i]]$p[j] * px[i]
        }
    }
    
    for(i in 1:m) {
        for(j in 1:n) {
            a <- (cpd[[i]]$p[j] * px[i])
            if(a > 0) {
                R <- R + a * (log2(cpd[[i]]$p[j]) - log2(qy[j]))
            }
        }
    }
    return(R);
}
