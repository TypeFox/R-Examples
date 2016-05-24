NCPUS <-
function(nchips = FALSE)
{
  
  if(nchips == FALSE) {
    nchips <- as.integer(detectCores()) 
    message("this platform has ", nchips, " cores.")
    return(NULL)
  }
  
  nchips <- as.integer(nchips)
  
  if(!is.integer(nchips))
    stop("'nchips' must be an integer.\n")
  
  if(nchips < 0)
    stop("'nchips' must be a positive integer if it is specified.\n")
  
  ans  <- makeCluster(nchips)
  return(ans)
}
