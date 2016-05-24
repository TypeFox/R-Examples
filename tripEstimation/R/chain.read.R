"chain.read" <-
function(filename) {
  con <- file(filename,open="rb")
  ## Read the dimensions of the array
  seek(con,0,origin="start",rw="r")
  dm <- readBin(con,"integer",3)
  ## Read the data and set dimensions
  n <- dm[1]*dm[2]*dm[3]
  A <- readBin(con,"double",n)
  if(length(A)!=n) stop("Read Error\n")
  close(con)
  array(A,dm)
}

