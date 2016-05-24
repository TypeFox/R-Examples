metri.p <- function(x1, x2, p = 2, w = 1){
           if (length(x1) != length(x2)) 
               stop("Dimensions of matrices do not match\n")
           if (length(w) == 1) 
               w = rep(w,length(x1))
           if (length(w) != length(x1)) 
               stop("Length of the vector does not match with the matrix\n")
           if (p > 0){
               (sum(w * (abs(x1 - x2) ^ p))) ^ (1/p)
           } 
           else {
               max(w * abs(x1 - x2))
           }
}
