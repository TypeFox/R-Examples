"unbound" <-
function(start, end=NULL) {
  if(is.matrix(start) && all(c("start", "end") %in% colnames(start))) {
     if(is.null(end)) end = start[, "end"]
     start = start[, "start"]
  }

  if(length(start)!=length(end))
    stop("start and end must are not the same length")
  ex <- NULL
  for(i in 1:length(start)) {
    ex <- c(ex, start[i]:end[i])
  }
  return(ex)
}

