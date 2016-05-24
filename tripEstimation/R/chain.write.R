`chain.write` <-
function(filename,A,append=FALSE) {
  A <- as.array(A)
  if(!append) {
    ## Open a binary connection for writing, truncating or
    ## creating the file if necessary
    con <- file(filename,open="wb")
    ## Write array dimensions
    writeBin(as.integer(dim(A)),con)
    ## Write the data and close
    writeBin(as.double(A),con)
    close(con)
  } else {
    ## Open a binary connection for appending
    con <- file(filename,open="r+b")
    ## Read the dimensions of the array
    seek(con,0,origin="start",rw="r")
    dm <- readBin(con,"integer",3)
    ## Check the data to append has compatible dimensions
    if(length(dm)!=length(dim(A)) || any(dm[-3]!=dim(A)[-3]))
      stop("Incompatible dimensions for appending\n")
    ## Append the new data
    seek(con,0,origin="end",rw="w")
    writeBin(as.double(A),con)
    ## Write new dimensions and close
    dm[3] <- dm[3]+dim(A)[3]
    seek(con,0,origin="start",rw="w")
    writeBin(as.integer(dm),con)
    close(con)
  }
}

