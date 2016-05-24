test.formatLong <- function() {
  df <- read.csv(system.file("csv","frog2001pcru.csv", package = "unmarked"))
  umf <- formatLong(df, type = "unmarkedFrameOccu")
  ## Add some assertions...
}



