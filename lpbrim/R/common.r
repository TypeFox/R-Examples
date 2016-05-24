#' @title Spread
#' @description Transforms values in v from m to M
#' @export
#' @param v Values to transform
#' @param m New min
#' @param M New max
spread = function (v, m = 0, M = 1)
{
   v <- v - min(v)
   v <- v/max(v)
   v <- v * (M - m)
   v <- v + m
   return(v)
}

#' @title Most frequent value
#' @description Returns the most frequent value in a vector
#' @export
#' @param vec The vector of labels
#' @param w The weights of each elements (currently not used)
mostFrequent = function(vec, w=NA)
{
   if(is.na(w)) w <- rep(1, length(vec)) # TODO use weights at some point?
   tvec <- table(vec)
   nvec <- as.vector(tvec[as.factor(vec)])
   return(which.max(nvec))
}

#' @title Edge list to adjacency matrix
#' @description From a three column matrix (from, to, value) to an adjacency matrix
#' @export
#' @param net Three-column matrix
netToAdjacency = function(net){
    x <- as.factor(net[,1])
    y <- as.factor(net[,2])
    z <- net[,3]
    adj <- matrix(0,
                  nrow=nlevels(x),
                  ncol=nlevels(y),
                  dimnames=list(levels(x), levels(y)))
    adj[cbind(x, y)] <- z
    return(adj)
}
