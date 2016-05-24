get.centroids <- function(map){
   if(! inherits(map,"bnd"))
      stop("Argument 'map' is not an object of class 'bnd'!")

   regions <- names(map)
   S <- length(regions)
  
   A <- vector(length=S)
   centroid <- matrix(nrow=S, ncol=2)

   for(k in 1:S){
      map.k <- map[[k]]
      nrows <- dim(map.k)[1]
      vec.A <- map.k[1:(nrows-1),1] * map.k[2:nrows,2] - map.k[2:nrows,1] * map.k[1:(nrows-1),2]
      A[k] <- 0.5 * sum(vec.A)
      vec.x <- (map.k[1:(nrows-1),1] + map.k[2:nrows,1]) * vec.A
      centroid[k,1] <- sum(vec.x)/(6*A[k])
      vec.y <- (map.k[1:(nrows-1),2] + map.k[2:nrows,2]) * vec.A
      centroid[k,2] <- sum(vec.y)/(6*A[k])
   }

   surrounding <- attr(map, "surrounding")
   whichAreInner <- which(sapply(surrounding, length) > 0L)
   for(l in seq_along(whichAreInner))
   {
       ## which is the inner polygon index?
       ind <- whichAreInner[l]
       
       ## which is the outer one?
       k <- which(regions == surrounding[[ind]])
       
       A[k] <- sign(A[k]) * (abs(A[k]) - abs(A[ind]))
         
       r <- c(centroid[ind,1] - centroid[k,1], centroid[ind,2] - centroid[k,2]) 
         a <- sqrt(r[1]^2 + r[2]^2)
         if(sum(r) != 0)
            r <- r/a
         s <- a*(-abs(A[ind])) / abs(A[k]) 

         centroid[k,] <- centroid[k,] + s * r
   }
   
   A <- abs(A)
   output <- cbind(A,centroid)
   rownames(output) <- regions
   colnames(output) <- c("area","centroid_x","centroid_y")
   return(output)
}
