print.epGraphs <-
function (x,...) {
	
  res.epGraphs <- x
  if (!inherits(res.epGraphs, "epGraphs")) stop ("no convenient data")
  cat("**ExPosition plotting data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(5, 2), list(1:5, c("name", "description")))
  
  res[1,] <- c("$fi.col","The colors for the row items.")
  res[2,] <- c("$fi.pch","The pch values for the row items.")  
  res[3,] <- c("$fj.col","The colors for the column items.")
  res[4,] <- c("$fj.pch","The pch values for the column items.")  
  res[5,] <- c("$constraints","Plotting constraints for axes.")  
  
  print(res)

}
