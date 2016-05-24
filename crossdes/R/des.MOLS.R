des.MOLS <- function(trt, k = trt) {
  
  primep100 <- c(2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 
    37, 41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97)
  if (!(trt %in% primep100)) {
    stop("trt is not a prime power less than 100.")
  }
  
  if (k%%1) {
    stop("k is not an integer.")
  } else {
    primefact <- matrix(c(2, 3, 2, 5, 7, 2, 3, 11, 13, 2, 17, 19, 23, 5, 3, 29, 
      31, 2, 37, 41, 43, 47, 7, 53, 59, 61, 2, 67, 71, 73, 79, 3, 83, 89, 97, 
      1, 1, 2, 1, 1, 3, 2, 1, 1, 4, 1, 1, 1, 2, 3, 1, 1, 5, 1, 1, 1, 1, 2, 
      1, 1, 1, 6, 1, 1, 1, 1, 4, 1, 1, 1), ncol = 2)
    trt <- primefact[primep100 == trt, ]
  }
  
  x <- MOLS(trt[1], trt[2])
  d <- NULL
  for (i in 1:(trt[1]^trt[2] - 1)) {
    d <- rbind(d, t(x[, , i]))
  }
  d[, 1:k]
} 
