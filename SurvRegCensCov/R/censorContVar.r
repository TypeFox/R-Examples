censorContVar <- function(x, LLOD = NA, ULOD = NA){
  
  n <- length(x)
  left <- rep(NA, n)
  right <- left
  
  if (is.na(LLOD) == FALSE){
    below <- which(x < LLOD)
    right[below] <- LLOD
  }
  
  if (is.na(ULOD) == FALSE){
    above <- which(x > ULOD)
    left[above] <- ULOD
  }
  
  bothmiss <- (apply(is.na(cbind(left, right)), 1, sum) == 2)
  left[bothmiss] <- x[bothmiss]
  right[bothmiss] <- x[bothmiss]
  
  res <- data.frame("left" = left, "right" = right)
  return(res)
}
