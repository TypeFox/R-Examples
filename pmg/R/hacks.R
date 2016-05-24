### These are also in gWidgetsRGtk/R/common.R -- in NAMESPACE NOW, along with untaintName,
#Paste = function(..., sep="", collapse="") {
#  x = unlist(list(...))
#  x = x[!is.na(x)]
#  x = x[x != "NA"]
#  paste(x, sep=sep, collapse=collapse)
#}

### ReadParseEvaL -- saves typing
#rpel = function(string, envir=.GlobalEnv) {
#  eval(parse(text=string), envir=envir)
#}

hack.as.data.frame = function(items) {
  ## check rectangular, or coerce to rectangules
  if(!(is.data.frame(items) || is.matrix(items) || is.vector(items))) {
    warning("Needs rectangular data, either a vector, matrix or data.frame")
    return(NA)
  }
  
  ## coerce to data frame
  if(is.vector(items)) {
    itemsName = deparse(substitute(items))
    items = data.frame(I(items))
    names(items) = itemsName
  }
  if(is.matrix(items)) {
    items = hack.as.data.frame.matrix(items) # fun in common.R
  }
  return(items)
}

## no easy way to not convert character vectors. This is a hack.
hack.as.data.frame.matrix = 
  function (x, row.names = NULL, optional = FALSE) 
  {
    d <- dim(x)
    nrows <- d[1]
    ir <- seq(length = nrows)
        ncols <- d[2]
    ic <- seq(length = ncols)
    dn <- dimnames(x)
    if (missing(row.names)) 
      row.names <- dn[[1]]
    collabs <- dn[[2]]
    if (any(empty <- nchar(collabs) == 0)) 
      collabs[empty] <- paste("V", ic, sep = "")[empty]
    value <- vector("list", ncols)
    if (mode(x) == "character") {
      for (i in ic) value[[i]] <- as.character(x[, i])
    }
    else {
      for (i in ic) value[[i]] <- as.vector(x[, i])
    }
    if (length(row.names) != nrows) 
      row.names <- if (optional) 
        character(nrows)
      else as.character(ir)
    if (length(collabs) == ncols) 
      names(value) <- collabs
    else if (!optional) 
      names(value) <- paste("V", ic, sep = "")
        attr(value, "row.names") <- row.names
    class(value) <- "data.frame"
    value
  }

