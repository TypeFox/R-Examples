print.mpGraphs <-
function (x ,...) {

  res.mpGraphs <- x
  if (!inherits(res.mpGraphs, "mpGraphs")) stop ("no convenient data")
  cat("**mexPosition plotting data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(3,2), list(1:3, c("name", "description")))
  
  res[1,] <- c('$table.col',"The colors of the tables")
  res[2,] <- c("$fi.col","The colors for the row items.")
  res[3,] <- c("$fj.col","The colors for the column items.")
  
  print(res)

}