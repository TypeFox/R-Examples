distClust <-
function(d, nobj.cluster, id.cluster)
{
   if(!inherits(d, "dist"))
      stop("'d' must be an object of class 'dist'!")
   dmat <- as.matrix(d)
   stopifnot(sum(nobj.cluster) == length(id.cluster))
   nc <- length(nobj.cluster)
   cl <- rep(1:nc, times = nobj.cluster)
   dc <- matrix(0, nc, nc)
   dimnames(dc) <- list(paste("cluster", 1:nc),
      paste("cluster", 1:nc))
   aux <- rbind(id.cluster, cl)
   diag(dmat) <- NA
   for(i in 1:nc) {
      for(j in 1:nc) {
         dc[i, j] <- mean(dmat[aux[1,][aux[2,] == i],
            aux[1,][aux[2,] == j]], na.rm = TRUE)
         if (is.nan(dc[i, j])) dc[i, j] <- 0
      }
   }
   return(dc)
}
