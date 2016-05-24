print.tinpoBoot <-
function (x,...) {

  res.tinpoBoot <- x
  if (!inherits(res.tinpoBoot, "tinpoBoot")) stop ("no convenient data")
  cat("**TInPosition Bootstrap output data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(2, 2), list(1:2, c("name", "description")))
  
  res[1,] <- c("$tests","Data for bootstrap ratio (BSR) tests.")
  res[2,] <- c("$boots","The array (rows, columns, depth) of bootstrap resampled & projected data.")
  
  print(res)

}
