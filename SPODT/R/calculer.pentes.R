calculer.pentes <-
function(data, data.sp)
{
    nbPts <- nrow(data.sp)
    pentes <- rep(0, nbPts*(nbPts-1)/2)
    ind1 <- rep(0, nbPts*(nbPts-1)/2)
    ind2 <- rep(0, nbPts*(nbPts-1)/2)
    nb.pentes.inf <- 0
    
    x <- data$x[match(data.sp$loc, data$loc)]
    y <- data$y[match(data.sp$loc, data$loc)]
    
    res <- .C("calculerPentes",
              as.integer(nbPts), as.integer(data.sp$loc), as.double(x), as.double(y),
              as.double(pentes), as.integer(ind1), as.integer(ind2), as.integer(nb.pentes.inf)
              
             )#ajouté ,PACKAGE="SPODT"
                      
    perm <- cbind.data.frame(res[[5]], res[[6]], res[[7]])
    colnames(perm) <- c("pente", "loc1", "loc2")
    perm <- rbind(perm)[order(perm$pente, perm$loc1, perm$loc2),]
    rownames(perm) <- 1:(nbPts*(nbPts-1)/2)

    if (res[[8]] != 0)
    {
        if (perm$pente[1] > 0)
        {
            perm$pente <- replace(perm$pente, (nbPts*(nbPts-1)/2-res[[8]]+1):(nbPts*(nbPts-1)/2), rep(perm$pente[1]/2, res[[8]]));
        }
        else
        {
            perm$pente <- replace(perm$pente, (nbPts*(nbPts-1)/2-res[[8]]+1):(nbPts*(nbPts-1)/2), rep(perm$pente[1]*2, res[[8]]));
        }
    }
     
    return(perm)
}
