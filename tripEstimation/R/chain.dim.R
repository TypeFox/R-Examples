`chain.dim` <-
function(filename) {
  con <- file(filename,open="rb")
  ## Read the dimensions of the array
  dm <- readBin(con,"integer",3)
  close(con)
  dm
}

