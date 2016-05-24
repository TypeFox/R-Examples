print.mexPosition.Output <- function (x ,...) {

  res.mexPosition.Output <- x
  if (!inherits(res.mexPosition.Output, "mexPosition.Output")) stop ("no convenient data")
  cat("**MexPosition output data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(2, 2), list(1:2, c("name", "description")))
  
  res[1,] <- c("$mexPosition.Data","All mexPosition classes output (data, factor scores, contributions, etc...)")
  res[2,] <- c("$Plotting.Data","All mexPosition & prettyGraphs plotting data (constraints, colors, etc...)")
  
  print(res)

}
