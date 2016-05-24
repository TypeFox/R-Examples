print.tinpoR2 <-
function (x,...) {

  res.tinpoR2 <- x
  if (!inherits(res.tinpoR2, "tinpoR2")) stop ("no convenient data")
  cat("**TInPosition R2 permutation test data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(3, 2), list(1:3, c("name", "description")))
  
  res[1,] <- c("$p.val","p-value of R2 permutation test.")
  res[2,] <- c("$r2.perm","The distribution of permuted R2 values.")
  res[3,] <- c("$r2","the observed r2 value.")  
  
  print(res)

}
