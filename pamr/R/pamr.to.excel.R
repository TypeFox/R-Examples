pamr.to.excel <- function(data, file, trace = TRUE) {
  if(is.null(data$x) | is.null(data$y) | is.null(data$genenames) | 
     is.null(data$geneid)) {
    stop("Invalid format for input data")
  }
  n <- nrow(data$x)
  p <- ncol(data$x)
  row1 <- paste("", "", sep = "\t")
  if(!is.null(data$samplelabels)) {
    for(j in 1:p) {
      row1 <- paste(row1, data$samplelabels[j], sep = "\t")
    }
    write(row1, file = file, append = FALSE)
  }
  row2 <- paste("", "", sep = "\t")
  if(!is.null(data$batchlabels)) {
    for(j in 1:p) {
      row2 <- paste(row2, data$batchlabels[j], sep = "\t")
    }
    write(row2, file = file, append = TRUE)
  }
  row3 <- paste("", "", sep = "\t")
  for(j in 1:p) {
    row3 <- paste(row3, data$y[j], sep = "\t")
  }
  write(row3, file = file, append = TRUE)
  for(i in 1:n) {
    if(trace) {
      cat(c("writing row number", i), fill = TRUE)
    }
    xx <- paste(data$gene.id[i], data$genenames[i],  sep = "\t")
    for(j in 1:ncol(data$x)) {
      xx <- paste(xx, data$x[i, j], sep = "\t")
    }
    write(xx, file = file, append = TRUE)
  }
  return()
}
