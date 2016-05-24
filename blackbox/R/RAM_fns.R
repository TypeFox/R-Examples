
get_os <- function() {
  if (.Platform$OS.type == "windows") { ## does not distinguish OSX (BSD) and mote standard unix-like systems
    "Windows"
  } else if (Sys.info()["sysname"] == "Darwin") { ## => OS X (or how to mess things)
    "OSX" 
  } else if (.Platform$OS.type == "unix") { 
    "unix"
  } else {
    stop("Unknown OS") ## but see http://conjugateprior.org/2015/06/identifying-the-os-from-r/
  }
}

RAMavail <- function() {
  ## windows case must be incorrect
  os <- get_os()
  if (os=="Windows") { 
    avail <- memory.size() ## FR->FR RAMused rather than available ?
  } else if (os=="OSX") { 
    bla <- system2("vm_stat",stdout=TRUE)[2] ## [1] => Free (but not Inactive)
    allpos <- gregexpr(" [0-9]",bla)[[1]] ## positions of blanks before numbers
    begpos <- 1L + allpos[length(allpos)]
    avail <- as.numeric(substr(bla,begpos,nchar(bla))) ## get last numeric value on the string
  } else { ## assuming unix-like system with a "free" function
    bla <- system2("free","-m",stdout=TRUE)[3] ## -m for value in Mo; [3] => "-/+ buffers/cache" line
    allpos <- gregexpr(" [0-9]",bla)[[1]] ## positions of blanks before numbers
    begpos <- 1L + allpos[length(allpos)]
    avail <- as.numeric(substr(bla,begpos,nchar(bla))) ## get last numeric value on the string
  }
  return(avail)
}

## function for memory management
#http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(all.names=TRUE,name = pos, pattern = pattern) ## FR 2015/09 name=pos, all=TRUE
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
# shorthands
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
#
lsOK <- function(which=1, n=10, ...) {
  ## still not fully effective, as e.g. lsOK(parent.frame()) will try (and fail) to evaluate ls.objects on each element of the evaluated frame...
  bla <- lapply(which,.ls.objects,..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
  if(length(which)>1L) {
    bla <- do.call(rbind,bla)
  } else bla <- bla[[1]]
  bla
}
