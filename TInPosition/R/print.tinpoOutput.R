print.tinpoOutput <-
function (x,...) {

  res.tinpoOutput <- x
  if (!inherits(res.tinpoOutput, "tinpoOutput")) stop ("no convenient data")
  cat("**TInPosition output data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(2, 2), list(1:2, c("name", "description")))
  
  res[1,] <- c("$Fixed.Data","All TExPosition descriptive data ($TExPosition.Data and $Plotting.Data).")
  res[2,] <- c("$Inference.Data","All TInPosition inference data (permutation, bootstrap tests).")
  
  print(res)

}
