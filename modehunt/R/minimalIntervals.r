minimalIntervals <- function(ints){

ints <- matrix(ints, ncol = 2)
nr <- length(ints[, 1])
res <- NULL

if (nr > 0){
     ints <- matrix(cbind(ints, 1:nr), ncol = 3)
     j.min <- NULL
     
     for (j in 1:nr){
         temp.j <- ints[j, ]
         subset <- (ints[, 1] >= temp.j[1]) * (ints[, 2] <= temp.j[2])
         
         if (length(subset) > 0){
             j.subsets <- matrix(ints[subset == 1, ], ncol = 3)
             j.length <- j.subsets[, 2] - j.subsets[, 1]
             if (length(j.length) > 0){j.min <- c(j.min, j.subsets[j.length == min(Inf, j.length), 3])}
         }
     }
     
     j.min <- sort(unique(j.min))     
     res <- ints[j.min, 1:2]
}

return(matrix(res, ncol = 2))
}
