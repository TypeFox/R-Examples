corRead <-
function(ref=NULL, names=NULL) {

  cat("\n")
  if (is.null(ref)) {
    ref <- file.choose()
    .dash(64)
    cat("File: \n")
    cat("   ", ref, "\n")
    .dash(64)
    cat("\n")
  }

  myc <- as.matrix(read.table(ref))
  if (!is.null(names)) colnames(myc) <- names
  rownames(myc) <- colnames(myc)

  return(myc)

}
