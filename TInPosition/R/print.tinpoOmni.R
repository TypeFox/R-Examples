print.tinpoOmni <-
function (x,...) {

  res.tinpoOmni <- x
  if (!inherits(res.tinpoOmni, "tinpoOmni")) stop ("no convenient data")
  cat("**TInPosition Omnibus Test data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(3, 2), list(1:3, c("name", "description")))
  
  res[1,] <- c("$p.val","p-value associated to omnibus permutation test.")
  res[2,] <- c("$inertia.perm","The distribution of permuted inertia values.")
  res[3,] <- c("$inertia","The observed inertia value.")  
  
  print(res)

}
