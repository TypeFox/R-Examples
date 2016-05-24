eD<-function(x,y){
   #Euclidean Distance
   #For multiple vectors in an array 
   #rows contain the k vectors
   #columns the n coordinates in the n-space
   #str(x) ==  matrix [1:k, 1:n]
   
   square <- (x-y)^2
   if(is.null(dim(square))){
        dim(square) <- c(1, length(square))
   }
   return(sqrt(rowSums(square)))
}

