print.tinpoComponents <-
function (x,...) {

  res.tinpoComponents <- x
  if (!inherits(res.tinpoComponents, "tinpoComponents")) stop ("no convenient data")
  cat("**TInPosition Components Permutation Test data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(3, 2), list(1:3, c("name", "description")))
  
  res[1,] <- c("$p.vals","p-values associated to components' permutation tests.")
  res[2,] <- c("$eigs.perm","The distributions of permuted eigenvalues.")
  res[3,] <- c("$eigs","The observed eigenvalues.")  
  
  print(res)

}
