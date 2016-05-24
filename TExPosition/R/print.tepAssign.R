print.tepAssign <-
function (x,...) {

  res.tepAssign <- x
  if (!inherits(res.tepAssign, "tepAssign")) stop ("no convenient data")
  cat("\n*The results of individual to group assignments are: \n\n")
  res <- array("", c(4, 2), list(1:4, c("name", "description")))
  
  res[1,] <- c("$r2","R-squared (between inertia/within inertia)")
  res[2,] <- c("$distances","Distances of observation to a group")
  res[3,] <- c("$assignment","Assignment of an observation to a group")
  res[4,] <- c("$confusion", "a confusion matrix (rows are estimates; columns are actual)")   
  
  print(res)

}
