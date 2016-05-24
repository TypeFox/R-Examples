pamr.from.excel <- function(file, ncols, sample.labels = FALSE, batch.labels = FALSE) {
  d <- scan(file, sep = "\t", what = "")
  dd <- matrix(d, ncol = ncols, byrow = TRUE)
  samplelabels <- NULL
  batchlabels <- NULL
  ii <- 1
  if(sample.labels) {
    samplelabels <- dd[1,  - (1:2)]
    ii <- ii + 1
  }
  if(batch.labels & !sample.labels) {
    batchlabels <- dd[1,  - (1:2)]
    ii <- ii + 1
  }
  if(batch.labels & sample.labels) {
    batchlabels <- dd[2,  - (1:2)]
    ii <- ii + 1
  }
  y <- dd[ii,  - (1:2)]
  geneid <- dd[ - (1:ii), 1]
  genenames <- dd[ - (1:ii), 2]
  x <- matrix(as.numeric(as.character(dd[ - (1:ii),  - (1:2)])), ncol = 
              ncols - 2)
  cat("",fill=TRUE)
  cat(c("Read in ", nrow(x), "genes"),fill=TRUE)
  cat(c("Read in ", ncol(x), "samples"),fill=TRUE)
  if(sample.labels){cat(c("Read in ", length(samplelabels), "sample labels"),fill=TRUE)}
  if(batch.labels){cat(c("Read in ", length(batchlabels), "batch labels"),fill=TRUE)}
  cat("",fill=TRUE)
  cat("Make sure these figures are correct!!", fill=TRUE)
  cat("",fill=TRUE)
  
  return(list(x = x, y = y, genenames = genenames, geneid = geneid, 
              samplelabels = samplelabels, batchlabels = batchlabels))
}

