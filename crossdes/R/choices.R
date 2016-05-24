choices <- function(trt, k = trt, maxsub = 1000) {
  
  choices <- logical(5)  # Which constructions work?
  
  b0 <- c((gamma(trt + 1)/gamma(trt - k + 1)), ifelse(!(trt%%2), trt, 2 * trt), 
    trt * (trt - 1), NA, trt * (trt - 1)/(k - 1))
  # How many subjects are required?
  
  primep100 <- c(2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 
    37, 41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97)  # Prime powers are not yet implemented. Get the first 35 ones manually.
  
  choices[1] <- b0[1] <= maxsub
  choices[2] <- (b0[2] <= maxsub) && (k == trt)
  choices[3] <- (b0[3] <= maxsub) && (trt %in% primep100[-1])
  
  
  patsub <- 0
  i <- 0
  ifelse(!(k%%2), maxi <- maxsub/k, maxi <- maxsub/(2 * k))
  while (!patsub && (i <= maxi)) {
    i <- i + 1
    if (!((i * k/trt)%%1) && !((i * k * (k - 1)/(trt * (trt - 1)))%%1)) {
      patsub <- (1 + (k%%2)) * i * k
    }
  }
  if ((k < trt) && (patsub > 0)) {
    choices[4] <- TRUE
    b0[4] <- patsub
  }
  
  
  if (((b0[5] <= maxsub) & (k < trt)) & ((((trt - 1)/(k - 1))%%1) == 0)) {
    choices[5] <- TRUE
  }
  
  out <- list(choices, b0)
  names(out) <- c("works", "number")
  out
} 
