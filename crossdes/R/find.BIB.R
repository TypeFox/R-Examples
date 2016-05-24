find.BIB <- function(trt, b, k, iter = 30) {
  i <- 0  # initialise stop criterions i, BIB
  BIB <- FALSE
  
  while ((!BIB) && (i < iter)) {
    i <- i + 1
    dummy <- optBlock(~., withinData = factor(1:trt), blocksizes = rep(k, b))
    des <- matrix(dummy$rows, byrow = TRUE, ncol = k)
    BIB <- all(isGYD(des, FALSE, FALSE)[[1]][1:4])
  }
  des
}

 
