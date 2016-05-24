#' @export
DStree.control <-
  function(minsplit = 20L, minbucket = round(minsplit/3), cp = 0.005,
           maxcompete = 4L,maxdepth=30,maxsurrogate=0)
  {
    maxsurrogate=0
    if (maxcompete < 0L) {
      warning("The value of 'maxcompete' supplied is < 0; the value 0 was used instead")
      maxcompete <- 0L
    }
   
    if (maxdepth > 30L) stop("Maximum depth is 30")
    if (maxdepth < 1L)  stop("Maximum depth must be at least 1")
    if (missing(minsplit) && !missing(minbucket)) minsplit <- minbucket * 3L
    
    
    
    list(minsplit = minsplit, minbucket = minbucket, cp = cp,
         maxcompete = maxcompete, maxsurrogate = maxsurrogate,
          maxdepth = maxdepth)
  }