myRound <- function(d){
 
    ind <- (d %% floor(d)) >= 0.5

    r <- floor(d)
    r[ind] <- ceiling(d[ind])

    return(r)}
