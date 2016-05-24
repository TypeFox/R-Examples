freq.calc <-
function(matr){    
     matr <- as.matrix(matr)
     if(!is.matrix(matr)){
         stop("The input data must be matrix!\n")
     }
     if(any(is.na(matr))){ 
         matr <- na.omit(matr)
         print(paste("NA found in matrix, and have been removed\n"))
     }
     matr[matr>1] <- 1
     result <- apply(matr, 2, sum)/nrow(matr)
     return(result)
}

