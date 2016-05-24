# This is file ../spam/R/foreign.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     




# Contains two sections:
# 1) Routines to transform spam objects to SparseM and Matrix
# 2) Functions to read (and write) MM and HB formats.



# 1a)  spam <-> SparseM
as.spam.matrix.csr <- function(x)
  {
#    if (is.matrix.csr(x)) {
      newx <- new("spam")
      slot(newx,"entries",check=FALSE) <- as.double( x@ra)
      slot(newx,"colindices",check=FALSE) <- x@ja
      slot(newx,"rowpointers",check=FALSE) <- x@ia
      slot(newx,"dimension",check=FALSE) <- x@dimension
      return(newx)
#    } else stop("Wrong object passed to 'as.spam.matrix.csr'")
  }


# The following should not be necessary because it is
# as."matrix.csr".spam and not "as.matrix".csr.spam.
# Is there anyway around this?
    
#as.matrix.csr.spam <- function(x,...) {
#  if (require('SparseM')){
#    newx <- new("matrix.csr")
#    slot(newx,"ra",check=FALSE) <- x@entries
#    slot(newx,"ja",check=FALSE) <- x@colindices
#    slot(newx,"ia",check=FALSE) <- x@rowpointers
#    slot(newx,"dimension",check=FALSE) <- x@dimension
#    return(newx)
#  }       
  
#}


# 1b) spam <-> Matrix

as.dgRMatrix.spam <- function(x) {
    if (requireNamespace('Matrix')) {
      newx <- new(p=0:0,'dgRMatrix')
      slot(newx,"x",check=FALSE) <- x@entries
      slot(newx,"j",check=FALSE) <- x@colindices-1L
      slot(newx,"p",check=FALSE) <- x@rowpointers-1L
      slot(newx,"Dim",check=FALSE) <- x@dimension
      return(newx)
    } 
  }

as.dgCMatrix.spam <- function(x)  {
    if (requireNamespace('Matrix')) {
      dimx <- x@dimension
      nz <- x@rowpointers[dimx[1] + 1] - 1
      z <- .Fortran("transpose", n = dimx[1], m = dimx[2],
                    a = as.double(x@entries),ja = x@colindices, ia = x@rowpointers,
                    entries = vector("double",nz), colindices = vector("integer", nz),
                    rowpointers = vector("integer", dimx[2] + 1),
                    NAOK = .Spam$NAOK,
                    PACKAGE = "spam")
      newx <- new(p=0:0,'dgCMatrix')
      slot(newx,"x",check=FALSE) <- z$entries
      slot(newx,"i",check=FALSE) <- z$colindices-1L
      slot(newx,"p",check=FALSE) <- z$rowpointers-1L
      slot(newx,"Dim",check=FALSE) <- dimx
      return(newx)
    } 
  }
    

as.spam.dgRMatrix <- function(x)  {
    
    if (is(x,'dgRMatrix')){
      if (identical(length(x@x),0L))  # zero matrix
        return(new("spam",rowpointers=c(1L,rep.int(2L,x@Dim[1])), dimension=x@Dim))

      newx <- new('spam')
      slot(newx,"entries",check=FALSE) <- x@x
      slot(newx,"colindices",check=FALSE) <- x@j+1L
      slot(newx,"rowpointers",check=FALSE) <- x@p+1L
      slot(newx,"dimension",check=FALSE) <- x@Dim
      return(newx)
    }
    stop("Wrong object passed to 'as.spam.dgRMatrix'")
  }
    
as.spam.dgCMatrix <- function(x)  {
    
    if (is(x,'dgCMatrix')){
      if (identical(length(x@x),0L))  # zero matrix
        return(new("spam",rowpointers=c(1L,rep.int(2L,x@Dim[1])), dimension=x@Dim))

      nz <- x@p[x@Dim[2] + 1]
      z <- .Fortran("transpose", n = x@Dim[2], m = x@Dim[1],
                    a = as.double(x@x),ja = x@i+1L, ia = x@p+1L,
                    entries = vector("double",nz), colindices = vector("integer", nz),
                    rowpointers = vector("integer", x@Dim[1] + 1),
                    NAOK = .Spam$NAOK,
                    PACKAGE = "spam")
      newx <- new('spam')
      slot(newx,"entries",check=FALSE) <- z$entries
      slot(newx,"colindices",check=FALSE) <- z$colindices
      slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
      slot(newx,"dimension",check=FALSE) <- x@Dim
      return(newx)
    }
    stop("Wrong object passed to 'as.spam.dgCMatrix'")
  }



# 2) Import and export
# taken from Matrix 0.999375-10 and adapted for spam
## Utilities for the Harwell-Boeing and MatrixMarket formats


readone <- function(ln, iwd, nper, conv)
# By Bates/Maechler from Matrix 0.999375-10 
{
    ln <- gsub("D", "E", ln)
    inds <- seq(0, by = iwd, length = nper + 1)
    (conv)(substring(ln, 1 + inds[-length(inds)], inds[-1]))
}

readmany <- function(conn, nlines, nvals, fmt, conv)
# By Bates/Maechler from Matrix 0.999375-10 
{
    if (!grep("[[:digit:]]+[DEFGI][[:digit:]]+", fmt))
	stop("Not a valid format")
    Iind <- regexpr('[DEFGI]', fmt)
    nper <- as.integer(substr(fmt, regexpr('[[:digit:]]+[DEFGI]', fmt), Iind - 1))
    iwd <- as.integer(substr(fmt, Iind + 1, regexpr('[\\.\\)]', fmt) - 1))
    rem <- nvals %% nper
    full <- nvals %/% nper
    ans <- vector("list", nvals %/% nper)
    for (i in seq_len(full))
	ans[[i]] <- readone(readLines(conn, 1, ok = FALSE),
			    iwd, nper, conv)
    if (!rem) return(unlist(ans))
    c(unlist(ans),
      readone(readLines(conn, 1, ok = FALSE), iwd, rem, conv))
}

read.HB <- function(file)
# Adapted from Bates/Maechler Matrix 0.999375-10 version
{
    if (is.character(file))
	file <- if (file == "") stdin() else file(file)
    if (!inherits(file, "connection"))
        stop("'file' must be a character string or connection")
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    hdr <- readLines(file, 4, ok = FALSE)
    Title <- sub('[[:space:]]+$', '', substr(hdr[1], 1, 72))
    Key <- sub('[[:space:]]+$', '', substr(hdr[1], 73, 80))
    totln <- as.integer(substr(hdr[2], 1, 14))
    ptrln <- as.integer(substr(hdr[2], 15, 28))
    indln <- as.integer(substr(hdr[2], 29, 42))
    valln <- as.integer(substr(hdr[2], 43, 56))
    rhsln <- as.integer(substr(hdr[2], 57, 70))
    if (!(t1 <- substr(hdr[3], 1, 1)) %in% c('C', 'R', 'P'))
        stop(paste("Invalid storage type:", t1))
    if (t1 != 'R') stop("Only numeric sparse matrices allowed")
    ## _FIXME: Patterns should also be allowed
    if (!(t2 <- substr(hdr[3], 2, 2)) %in% c('H', 'R', 'S', 'U', 'Z'))
        stop(paste("Invalid storage format:", t2))
    if (!(t3 <- substr(hdr[3], 3, 3)) %in% c('A', 'E'))
        stop(paste("Invalid assembled indicator:", t3))
    nr <- as.integer(substr(hdr[3], 15, 28))
    nc <- as.integer(substr(hdr[3], 29, 42))
    nz <- as.integer(substr(hdr[3], 43, 56))
    nel <- as.integer(substr(hdr[3], 57, 70))
    ptrfmt <- toupper(sub('[[:space:]]+$', '', substr(hdr[4], 1, 16)))
    indfmt <- toupper(sub('[[:space:]]+$', '', substr(hdr[4], 17, 32)))
    valfmt <- toupper(sub('[[:space:]]+$', '', substr(hdr[4], 33, 52)))
    rhsfmt <- toupper(sub('[[:space:]]+$', '', substr(hdr[4], 53, 72)))
    if (!is.na(rhsln) && rhsln > 0) {
        h5 <- readLines(file, 1, ok = FALSE)
    }
    ptr <- readmany(file, ptrln, nc + 1, ptrfmt, as.integer)
    ind <- readmany(file, indln, nz, indfmt, as.integer)
    vals <- readmany(file, valln, nz, valfmt, as.numeric)

    # Spam related changes:
    if (t3 =="E")
        stop("Only assembled Harwell-Boeing formats implemented")      
    z <- .Fortran("transpose", n = nc, m = nr,
                  a = vals,ja = ind, ia = ptr,
                  entries = vector("double",nz), colindices = vector("integer", nz),
                  rowpointers = vector("integer", nr + 1),
                  NAOK = .Spam$NAOK,
                  PACKAGE = "spam")
    newx <- new('spam')
    slot(newx,"entries",check=FALSE) <- z$entries
    slot(newx,"colindices",check=FALSE) <- z$colindices
    slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
    slot(newx,"dimension",check=FALSE) <- c(nr, nc)
    if (t2 %in% c('H', 'S'))
      newx <- newx+t.spam(newx)-diag.spam(spam(newx))
    if (t2 =="Z")
      newx <- newx-t.spam(newx)
    return(newx)
}



# alternatives are implementing
# http://math.nist.gov/MatrixMarket/mmio/f/mmiof77.html
read.MM <- function(file)  {
  if (is.character(file))
    file <- if(file == "") stdin() else file(file)
  if (!inherits(file, "connection"))
    stop("'file' must be a character string or connection")
  if (!isOpen(file)) {
    open(file)
    on.exit(close(file))
  }
  scan1 <- function(what, ...)
    scan(file, nmax = 1, what = what, quiet = TRUE, ...)
  
  if ((hdr <- tolower(scan1(character()))) != "%%matrixmarket")  # RF: added a to lower
    stop("file is not a MatrixMarket file")
  if (!(typ <- tolower(scan1(character()))) %in% "matrix")
    stop("type '", typ, "' not recognized")
  if (!(repr <- tolower(scan1(character()))) %in% c("coordinate", "array"))
    stop("representation '", repr, "' not recognized")
  elt <- tolower(scan1(character()))
  if (!elt %in% c("real", "complex", "integer", "pattern"))
    stop("element type '", elt, "' not recognized")
  sym <- tolower(scan1(character()))
  if (!sym %in% c("general", "symmetric", "skew-symmetric", "hermitian"))
    stop("symmetry form '", sym, "' not recognized")
  nr <- scan1(integer(), comment.char = "%")
  nc <- scan1(integer())
    # code from now on differs from Matrix one...
  if (repr == "coordinate") {
    nz <- scan1(integer())
    switch(elt,
           "real"    = { what <- list(i= integer(), j= integer(), x= numeric())},
           "integer" = { what <- list(i= integer(), j= integer(), x= numeric())
                         warning("'integer' format coerced to 'double'", call. = FALSE)         },
           "pattern" = { what <-  list(i= integer(), j= integer())
                         warning("matrix elements assumed as 1 ('pattern' format)", call. = FALSE)    },
           "complex" = { what <- list(i= integer(), j= integer(), x= numeric(), y= numeric())
                         warning("retaining only real part of 'complex' format", call. = FALSE) }           )
    
    z <- scan(file, nmax = nz, quiet = TRUE, what= what)
    newx <- spam.list(list(ind=cbind(z$i,z$j),x= if(elt=="pattern") rep.int(1,nz) else z$x ), nr,nc)
    
    if (sym %in% c("symmetric", "hermitian"))  {
      dim(newx) <- rep(max(nr,nc),2)
      newx <- newx+t.spam(newx)-diag.spam(diag(newx))
    }
    if (sym=="skew-symmetric") {
      dim(newx) <- rep(max(nr,nc),2)
      newx <- newx-t.spam(newx)
    }
  }  else {
    nz <- nr*nc
    x <- scan(file, nmax = nz, quiet = TRUE, what=numeric())
    z <- .Fortran("spamdnscsr", nrow = nr, ncol = nc,
                  x = x, nr, entries = vector("double",nz),
                  colindices = vector("integer", nz), rowpointers = vector("integer",nr + 1),
                  eps = spam.options('eps'), NAOK = TRUE,
                  PACKAGE = "spam")
    
    warning("returning a (possibly) dense 'spam' object", call. = FALSE)
    nz <- z$rowpointers[nr+1]-1
    if (identical(nz, 0L))
      return(new("spam",rowpointers=c(1L,rep.int(2L,nr)), dimension=c(nr,nc)))
  
    newx <- new("spam")
    slot(newx,"entries",check=FALSE) <- z$entries[1:nz]
    slot(newx,"colindices",check=FALSE) <- z$colindices[1:nz]
    slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
    slot(newx,"dimension",check=FALSE) <- c(nr,nc)
  }
  return(newx)
}
