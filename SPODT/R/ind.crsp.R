ind.crsp <-
function(crsp, loc1, loc2)
{
    nLoc <- length(crsp)
    ind1 <- rep(0, length(loc1))
    ind2 <- rep(0, length(loc2))
    res <- .C("indCrsp",
              as.integer(nLoc), as.integer(crsp), as.integer(loc1), as.integer(loc2),
              as.integer(ind1), as.integer(ind2)
              
             )#ajout ,PACKAGE="SPODT"

    return(list(i1=res[[5]], i2=res[[6]]))
}
