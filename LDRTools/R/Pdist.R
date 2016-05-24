Pdist <- function(x, weights = "constant")
    {
    weights <- match.arg(weights, c("constant", "inverse", "sq.inverse"))
    
    n <- length(x)
    k <- sapply(x, matrix.trace)
    
    weight <- switch(weights, constant = {ifelse(k==0, 0, 1)}
                            , inverse = {ifelse(k==0, 0, 1/k)}
                            , sq.inverse = {ifelse(k==0, 0, 1/sqrt(k))}
                            )
                            
    x.weighted <-  mapply("*", x, weight, SIMPLIFY=FALSE)
    
    Pdist <- numeric(n*(n-1)/2)
    ii <- 0
    for (i in (2:n))
        {
        for  (j in (1:(i-1)))
            {
            ii <- ii+1
            Pdist[ii] <- norm((x.weighted[[i]]-x.weighted[[j]]), type="F")^2
            }
        }
    Pdist <- Pdist/2
    
    attr(Pdist, "Size") <- n
    attr(Pdist, "Labels") <- names(x)
    attr(Pdist, "Diag") <- FALSE
    attr(Pdist, "Upper") <- FALSE
    attr(Pdist, "methods") <- weights
    
    class(Pdist) <- "dist"
    Pdist
    }
