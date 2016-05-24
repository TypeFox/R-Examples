# S3 accessor methods
extractFUN <- function(x) {
  UseMethod("extractFUN")
}
`extractFUN<-` <- function(x, value) {
  UseMethod("extractFUN<-")
}
extractFUN.mmap <- function(x) {
  x$extractFUN
}
`extractFUN<-.mmap` <- function(x, value) {
  x$extractFUN <- value
  x
}
replaceFUN <- function(x) {
  UseMethod("replaceFUN")
}
`replaceFUN<-` <- function(x, value) {
  UseMethod("replaceFUN<-")
}
replaceFUN.mmap <- function(x) {
  x$replaceFUN
}
`replaceFUN<-.mmap` <- function(x, value) {
  x$replaceFUN <- value
  x
}

# Basic S3 methods
head.mmap <- function(x, n=6L, ...) {
  x[1:(min(NROW(x),n))]
}
tail.mmap <- function(x, n=6L, ...) {
  x[(NROW(x)-n):NROW(x)]
}

str.mmap <- function (object, ...) 
{
    print(object)
    cat("  data         :") 
    str(object$data)
    cat("  bytes        :")
    str(object$bytes)
    cat("  filedesc     :")
    str(object$filedesc)
    cat("  storage.mode :") 
    str(object$storage.mode)
    cat("  pagesize     :")
    str(object$pagesize)
    cat("  dim          :")
    print(object$dim)
}

summary.mmap <- function(object) { str(object) }

print.mmap <- function(x, ...) {
  stopifnot(is.mmap(x))
  file_name <- names(x$filedesc)
  if(nchar(file_name) > 10)
    file_name <- paste(substring(file_name,0,10),"...",sep="")
  type_name <- switch(typeof(x$storage.mode),
                      "list"="struct",
                      "integer"="int",
                      "logical"="logi",
                      "double"="num",
                      "complex"="cplx",
                      "character"="chr",
                      "raw"="raw")
  if(type_name == "struct") {
    firstN <- x[1][[1]]
  } else {
  firstN <- head(x)
  firstN <- if(cumsum(nchar(firstN))[length(firstN)] > 20) {
                firstN[1:min(3,length(x))]
              } else {
                firstN
              }
  }
  if( !is.null(x$dim)) { # has dim
  cat(paste("<mmap:",file_name,">  (",class(x$storage.mode)[2],") ",
            type_name," [1:", nrow(x),", 1:", ncol(x),"]",sep=""),firstN,"...\n")
  } else {
  cat(paste("<mmap:",file_name,">  (",class(x$storage.mode)[2],") ",
            type_name," [1:", length(x),"]",sep=""),firstN,"...\n")
  }
}
print.summary_mmap <- function() {}

close.mmap <- function(con, ...) {
  munmap(con)
}

# creat flags using upper case symbols/strings
# mmapFlags(PROT_READ,PROT_WRITE) OR mmapFlags(PROT_READ | PROT_WRITE)
mmapFlags <- function(...) {
  flags <- as.character(match.call(call=sys.call())[-1])
  if(nargs()==1) {
    flags <- gsub(" ","",unlist(strsplit(flags,"\\|")))
    flags <- gsub('\"',"",flags) # in case "" | ""
  }
  .Call("mmap_mkFlags", flags, PKG="mmap")
}

# S3 constructor
mmap <- function(file, mode=int32(), 
                 extractFUN=NULL, replaceFUN=NULL,
                 prot=mmapFlags("PROT_READ","PROT_WRITE"),
                 flags=mmapFlags("MAP_SHARED"),len,off=0L,
                 ...) {
    if(missing(file))
      stop("'file' must be specified")
    # pageoff is the offset from page boundary
    # off is the page-aligned offset
    #   e.g. off=22 would set off=0 and pageoff=22 on a system with 4096 page sizing
    pageoff <- off %% pagesize()
    off <- off - pageoff
    if(missing(len))
      len <- file.info(file)$size - off - pageoff
    
    mmap_obj <- .Call("mmap_mmap", 
                      as.Ctype(mode),
                      file,
                      as.integer(prot), 
                      as.integer(flags), 
                      as.double(len),
                      as.integer(off),
                      as.integer(pageoff),
                      PKG="mmap")
    reg.finalizer(mmap_obj, mmap_finalizer, TRUE)
    mmap_obj$filedesc <- structure(mmap_obj$filedesc, .Names=file)
    mmap_obj$extractFUN <- extractFUN
    mmap_obj$replaceFUN <- replaceFUN
    class(mmap_obj) <- "mmap"
    return(mmap_obj)
}

# S3 destructor
munmap <- function(x) {
  if(!is.mmap(x))
    stop("mmap object required to munmap")
  invisible(.Call("mmap_munmap", x, PKG="mmap"))
}

mmap_finalizer <- function(x) {
  if(is.mmap(x))
    invisible(.Call("mmap_munmap", x, PKG="mmap"))
}

msync <- function(x, flags=mmapFlags("MS_ASYNC")) {
  if(!is.mmap(x))
    stop("mmap object required to munmap")
  .Call("mmap_msync", x, as.integer(flags), PKG="mmap")
}

mprotect <- function(x, i, prot) {
  # i indicates the start and length of protection

  # TODO: add ability to protect multiple pages in a
  # range
  .Call("mmap_mprotect", x, i, prot, PKG="mmap")
}

is.mmap <- function(x) {
  inherits(x, "mmap") && .Call("mmap_is_mmapped",x,PKG="mmap")
}

`[.mmap` <- function(x, i, j, ...) {
  if(!x$bytes) stop('no data to extract')
  if( is.struct(x$storage.mode) || is.null(x$dim)) {
    if(missing(i))
      i <- 1:length(x)
    if(missing(j))
      j <- 1:length(x$storage.mode)
    if(is.character(j))
      j <- match(j, names(x$storage.mode))
    DIM <- NULL
  } else {
    if( missing(i))
      i <- 1:dim(x)[1]
    if( missing(j))
      j <- 1:dim(x)[2]
    DIM <- c(length(i),length(j))
    i <- .Call("convert_ij_to_i", as.double(dim(x)[1]), as.integer(i), as.integer(j))
    j <- 1L
  }
  j <- j[j>0] # only positive values
  xx <- .Call("mmap_extract", i, as.integer(j), DIM, x, PKG="mmap")
  names(xx) <- names(x$storage.mode)[j]
  if(is.null(extractFUN(x))) {
    xx
  } else as.function(extractFUN(x))(xx)
}

`[<-.mmap` <- function(x, i, j, ..., sync=TRUE, value) {
  # add type checking/coercing at the C-level
  if(!x$bytes) stop('no data to extract')
  if(is.struct(x$storage.mode) && !is.list(value))
    value <- list(value)
  if( is.struct(x$storage.mode) || is.null(x$dim)) {
    if(missing(i))
      i <- 1:length(x)
    if(missing(j))
      j <- 1:length(x$storage.mode)
    if(is.character(j))
      j <- match(j, names(x$storage.mode))
    if(length(i) != length(value))
      if(is.list(value))
        value <- lapply(value, rep, length.out=length(i))
      else value <- rep(value, length.out=length(i))
  } else { # has dimension
    if(missing(i))
      i <- 1:dim(x)[1]
    if(missing(j))
      j <- 1:dim(x)[2]
    if(is.character(j))
      j <- match(j, names(x$dimnames))
    i <- .Call("convert_ij_to_i", as.double(dim(x)[1]), as.integer(i), as.integer(j))
    j <- 1L
    if(length(i) != length(value))
      value <- rep(value, length.out=length(i))
  }
# likely we need to check for list()/struct to correctly handle in C
  if(max(i) > length(x) || min(i) < 0)
    stop("improper 'i' range")
  .Call("mmap_replace", i, j, value, x, PKG="mmap") 
  if(sync)
    msync(x)
  x
}

length.mmap <- function(x) {
  size_in_bytes <- x$bytes
  size <- attr(x$storage.mode,"bytes")
  if( class(x$storage.mode)[2] == 'bits')
    trunc(size_in_bytes/size) * 32L
  else
  trunc(size_in_bytes/size)
}

`length<-.mmap` <- function(x, value) {
  # should have some mechanism to sanity check length
  # as to not pass end of file
  size_in_bytes <- value * attr(x$storage.mode, "bytes")
  if(file.info(names(x$filedesc))$size < size_in_bytes) {
    stop("cannot increase an mmap file's size") # implement something automatic here
  }
  x$bytes <- as.double(size_in_bytes)
  x 
}


# coerce to disk object and mmap back in.  Need
# to register a proper finalizer in the C code, else
# we are likely to end up with memory leaks.  For now
# this is not too probable, and not too dangerous.
# Famous last words ...

as.mmap <- function(x, mode, file,...) {
  UseMethod("as.mmap")
}

#as.mmap.data.frame <- function(x, mode=as.struct(x), file, ...) 

as.mmap.raw <- function(x, mode=raw(), file=tempmmap(), ...) {
  writeBin(x, file)
  mmap(file, raw())
}

as.mmap.integer <- function(x,
                            mode=integer(),
                            file=tempmmap(),
                            ...) {
  nbytes <- attr(as.Ctype(mode),"bytes")
  if(nbytes == 3) {
    writeBin(writeBin(x,raw())[1:(length(x)*4) %% 4 != 0], file)
  } else writeBin(x, file, size=nbytes)
  mmap(file, as.Ctype(mode))
}
as.mmap.double <- function(x,
                            mode=double(),
                            file=tempmmap(),
                            ...) {
  nbytes <- attr(as.Ctype(mode),"bytes")
  writeBin(x, file, size=nbytes)
  mmap(file, as.Ctype(mode))
}
as.mmap.complex <- function(x,
                            mode=complex(),
                            file=tempmmap(),
                            ...) {
  nbytes <- attr(as.Ctype(mode),"bytes")
  writeBin(x, file, size=nbytes)
  mmap(file, as.Ctype(mode))
}

as.mmap.character <- function(x, 
                              mode=char(nchar(x[1])), 
                              file=tempmmap(), force=FALSE, ...) {
  if( !all(nchar(x) == nchar(x[1]))) {
    if(!force)
      stop("requires fixed-width character vector. Use make.fixedwidth first.")
    x <- make.fixedwidth(x)
  }
  #if( !identical(mode, char(nchar(x[1])))){
  writeBin(x, file, size = nbytes)
  mmap(file, as.Ctype(mode))
}

as.mmap.matrix <- function(x, mode, file=tempmmap(), ...) {
  if( missing(mode))
    mode <- vector(storage.mode(x))
  DIM <- dim(x)
  dim(x) <- NULL
  x <- as.mmap(x, mode, file, ...)
  dim(x) <- DIM
  x
}

as.mmap.data.frame <- function(x, mode, file, ...) {
  if( !missing(mode))
    warning("'mode' argument currently unsupported")
  tmp <- tempfile()
  bytes.needed <- function(x) {
    sum(sapply(x, function(X) sizeof(as.Ctype(do.call(storage.mode(X),list(0)))))) * NROW(x)
  }
  writeBin(rep(as.raw(0), bytes.needed(x)), tmp)
  struct.type <- do.call(struct, sapply(x, function(X) (do.call(storage.mode(X),list(0)))))
  m <- mmap(tmp, struct.type, extractFUN=function(x) as.data.frame(x))
  dimnames(m) <- list(NULL, colnames(x))
  for(i in 1:NCOL(x)) {
    m[,i] <- x[,i]
  }
  m
}


tempmmap <- function(tmpdir=tempdir()) {
  tempfile("mmap",tmpdir)
}

pagesize <- function() {
  .Call("mmap_pagesize")
}
