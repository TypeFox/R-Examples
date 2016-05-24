# CHNOSZ/util.character.R
# functions to work with character objects

c2s <- function(x, sep=' ') {
  # make a string out of a character vector
  if(length(x) %in% c(0,1)) return(x)
  s <- paste(x,collapse=sep)
  return(s)
}

s2c <- function(x,sep=NULL,keep.sep=TRUE) {
  # recursively split 'x' according to separation strings in 'sep'
  do.split <- function(x,sep,keep.sep=TRUE) {
    # split the elements of x according to sep
    # output is a list the length of x
    if(is.list(x)) stop("x is a list; it must be a character object (can have length > 1)")
    x <- as.list(x)
    for(i in 1:length(x)) {
      # do the splitting
      xi <- strsplit(x[[i]],sep,fixed=TRUE)[[1]]
      # paste the separation term term back in
      if(keep.sep & !is.null(sep)) {
        xhead <- character()
        xtail <- xi
        if(length(xi) > 1) {
          xhead <- head(xi,1)
          xtail <- tail(xi,-1)
          # in-between matches
          xtail <- paste("",xtail,sep=sep)
        } 
        # a match at the end ... grep here causes problems
        # when sep contains control characters (e.g. protein.refseq)
        #if(length(grep(paste(sep,"$",sep=""),x[[i]]) > 0)) xtail <- c(xtail,sep)
        # use substr instead
        nx <- nchar(x[[i]])
        ns <- nchar(sep)
        if(substr(x[[i]],nx-ns+1,nx) == sep) xtail <- c(xtail,sep)
        xi <- c(xhead,xtail)
      }
      x[[i]] <- xi
    }
    return(x)
  }
  # now do it!
  for(i in 1:length(sep)) x <- unlist(do.split(x,sep[i],keep.sep=keep.sep))
  return(x)
}

can.be.numeric <- function(x) {
  # return FALSE if length of argument is zero
  if(length(x)==0) return(FALSE)
  if(length(x)>1) return(as.logical(sapply(x,can.be.numeric)))
  # don't warn about NAs in as.numeric
  oldopt <- options(warn=-1)
  cb <- FALSE
  if(!is.na(as.numeric(x))) cb <- TRUE
  if(x %in% c('.','+','-')) cb <- TRUE
  # let the user have their way
  options(oldopt)
  return(cb)
}

