get.segments <- function(i1, i2 = NULL, oob.size = 1, max.seg = 100)
{
  if (is.null(i2)) {
    Y <- i1
    
    i1 <- which(Y == names(table(Y))[1])
    i2 <- which(Y == names(table(Y))[2])
  }
  
  n1 <- length(i1)
  n2 <- length(i2)

  if (oob.size == 1) {
    segments <- rbind(rep.int(i1, n2), rep.int(i2, rep(n1, n2)))
    if (!is.null(max.seg) & ncol(segments) > max.seg) 
      segments <- segments[,sample(ncol(segments), max.seg)]
  } else {
    if (is.null(max.seg))
      stop("max.seg cannot be NULL when oob.size is larger than 1")
    
    segments <- sapply(1:max.seg,
                       function(i) c(sample(i1, oob.size),
                                     sample(i2, oob.size)))
  }

  segments
}
