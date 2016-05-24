read.yearmon <- function(file, format = "", tz = "", FUN = NULL,
                         regular = FALSE, index.column = 1, ...)
{
  ## `file' and `...' are simply passed to read.table
  ## the first column is interpreted to be the index, the rest the coredata
  ## it is transformed to an arbitrary index class by `FUN'
  ## defaults for `FUN' are guessed and are numeric, Date or POSIXct

  ## read data
  rval <- read.table(file, ...)

  ## if `file' does not contain data
  Nr <- NROW(rval)
  if(Nr < 1) {
    if(is.data.frame(rval)) rval <- as.matrix(rval)
    if(NCOL(rval) > 1) rval <- rval[,-index.column]
    rval <- zoo(rval)
    return(rval)
  }

  ## extract index
  if(NCOL(rval) < 1) stop("data file must specify at least one column")
  
  ## extract index, retain rest of the data
  if (NCOL(rval) == 1) ix <- seq(length = Nr)
  else {
    ix <- rval[,index.column]
    rval <- rval[,-index.column]
  }
  if(is.factor(ix)) ix <- as.character(ix)
  if(is.data.frame(rval)) rval <- as.matrix(rval)
    
  ## index transformation functions
  toDate <- if(missing(format)) function(x) as.Date(as.character(x))
              else function(x) as.Date(as.character(x), format = format)
  toPOSIXct <- function(x) as.POSIXct(as.character(x), tz = tz)
  toDefault <- function(x) {
    rval <- try(toPOSIXct(x), silent = TRUE)
    if(inherits(rval, "try-error"))
      rval <- try(toDate(x), silent = TRUE)
    else {
      hms <- as.POSIXlt(rval)
      hms <- hms$sec + 60 * hms$min + 3600 * hms$hour
      if(isTRUE(all.equal(hms, rep.int(hms[1], length(hms))))) {
        rval2 <- try(toDate(x), silent = TRUE)
        if(!inherits(rval2, "try-error")) rval <- rval2
      }
    }
    if(inherits(rval, "try-error")) rval <- rep(NA, length(x))
    return(rval)
  }
  toNumeric <- function(x) x
  
  ## setup default FUN
  if(is.null(FUN)) {
    FUN <- if(!missing(format)) toDate
           else if(!missing(tz)) toPOSIXct
           else if(is.numeric(ix)) toNumeric
           else toDefault        
  }
  
  ## compute index from (former) first column
  ix <- FUN(ix)
  
  ## sanity checking
  if(any(is.na(ix))) stop("index contains NAs")
  if(length(ix) != Nr) stop("index does not match data")

  ## Can convert ix to yearmon?
  ymx <- as.yearmon(ix)
  {
    nMo <- length(unique(ymx))
    if(nMo == Nr) 
      names(ymx) <- ix
    else{
      warning("'yearmon' index not allowed:  Only ", nMo,
              " different months in ", Nr, " records read. ",
              " Returning a zoo object with the dates as read.")
      ymx <- ix
    }
  }
  ## setup zoo object and return
  rval <- zoo(rval, ymx)
  if(regular && is.regular(rval)) rval <- as.zooreg(rval)
  return(rval)
}
