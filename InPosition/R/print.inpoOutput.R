print.inpoOutput <-
function (x,...) {

  res.inpoOutput <- x
  if (!inherits(res.inpoOutput, "inpoOutput")) stop ("no convenient data")
  cat("**InPosition output data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(2, 2), list(1:2, c("name", "description")))
  
  res[1,] <- c("$Fixed.Data","All ExPosition descriptive data ($ExPosition.Data and $Plotting.Data).")
  res[2,] <- c("$Inference.Data","All InPosition inference data (permutation, bootstrap tests).")
  
  print(res)

}
