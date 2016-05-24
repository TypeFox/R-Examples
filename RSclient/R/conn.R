RSconnect <- function(host="localhost", port=6311) {
  c <- socketConnection(host,port,open="a+b",blocking=TRUE)
  a <- readBin(c,"raw",32)
  if (!length(a)) { close(c); stop("Attempt to connect to Rserve timed out, connection closed") }
  if (length(a) != 32 || !length(grep("^Rsrv01..QAP1",rawToChar(a))))
    stop("Invalid response from Rserve")
  return( c )
}

RSeval <- function(c, expr) {
  r <- if (is.character(expr)) serialize(parse(text=paste("{",paste(expr,collapse="\n"),"}"))[[1]],NULL) else serialize(expr, NULL)
  writeBin(c(0xf5L, length(r), 0L, 0L), c, endian="little")
  writeBin(r, c)
  b <- readBin(c,"int",4,endian="little")
  if (length(b)<4 || b[1] != 65537L) stop("remote evaluation failed")
  unserialize(readBin(c,"raw",b[2]))
}

RSassign <- function (c, obj, name = deparse(substitute(obj))) {
  r <- serialize(list(name, obj), NULL)
  writeBin(c(0xf6L,length(r),0L,0L), c, endian="little")
  writeBin(r, c)
  b <- readBin(c,"int",4,endian="little")
  if (length(b)<4 || b[1] != 65537L)
    stop("remote assign failed")
  invisible(obj)
}

RSclose <- function(c) close(c)

# convert an array of unsigned integers into raw verctor safely
# by converting 16-bits at a time
.safe.int <- function(data) {
  r <- raw(length(data) * 4)
  j <- 1
  for (i in data) {
    hi <- as.integer(i / 0x10000 + 0.5)
    lo <- as.integer( (i - hi*0x10000) + 0.5)
    rs <- writeBin(c(lo, hi), raw(), endian="little")
    r[j] <- rs[1]
    r[j+1] <- rs[2]
    r[j+2] <- rs[5]
    r[j+3] <- rs[6]
    j <- j + 4
  }
  r
}

RSdetach <- function( c ) RSevalDetach( c, "" )

RSevalDetach <- function( c, cmd="" ) {
  # retrieve the host name from the connection (possibly unsafe!)
  host <- substr(strsplit(summary(c)$description,":")[[1]][1],3,999)
  if ( cmd != "" ) {
    r <- paste("serialize({", cmd[1], "},NULL)")
    l <- nchar(r[1])+1
    writeBin(as.integer(c(0x031,l+4,0,0,4+l*256)), c, endian="little")
    writeBin(as.character(r[1]), c)
    b <- readBin(c,"int",4,endian="little")
    if (b[1]%%256 == 2 || b[2] < 12) stop("Eval/detach failed with error: ",b[1]%/%0x1000000)
    ## We don't need "isLarge" because we never get large data back
  } else {
    l <- 0
    writeBin(as.integer(c(0x030,l+4,0,0,4+l*256)), c, endian="little")
    b <- readBin(c,"int",4,endian="little")
    if (b[1]%%256 != 1) stop("Detach failed with error: ",b[1]%/%0x1000000)
  }
  msgLen <- b[1]%/%256
  a <- readBin(c,"int",2,signed=FALSE,endian="little")
  if (!length(a)) { close(c); stop("Rserve connection timed out and closed") }
  ## a[1] is DT_INT, a[2] is the payload (port#)
  port <- a[ 2 ]
  readBin(c,"raw",4) ## this should be DT_BYTESTREAM
  key <- readBin(c,"raw",msgLen-12)
  RSclose(c)
  list( port=port, key=key, host=host )
}

RSattach <- function(session) {
  c <- socketConnection(session$host,session$port,open="a+b",blocking=TRUE)
  writeBin( session$key, c )
  b <- readBin(c,"int",4,endian="little")
  if (!length(b)) { close(c); stop("Rserve connection timed out and closed") }
  if (b[1]%%256 != 1) stop("Attach failed with error: ",b[1]%/%0x1000000)
  c
}

RSlogin <- function(c, user, pwd, silent=FALSE) {
  r <- paste(user,pwd,sep="\n")
  l <- nchar(r[1])+1
  writeBin(as.integer(c(1,l+4,0,0,4+l*256)), c, endian="little")
  writeBin(as.character(r[1]), c)
  b <- readBin(c,"int",4,endian="little")
  if (!length(b)) { close(c); stop("Rserve connection timed out and closed") }
  ##cat("header: ",b[1],", ",b[2],"\n")    
  msgLen <- b[2]
  if (msgLen > 0) a <- readBin(c,"raw",msgLen)
  if (b[1]%%256 != 1 && !silent) stop("Login failed with error: ",b[1]%/%0x1000000)
  invisible(b[1]%%256 == 1)
}

RSserverEval <- function(c, expr) {
  if (is.language(expr)) expr <- deparse(expr)
  if (!is.character(expr)) stop("expr must me a character vector, name, call or an expression")
  r <- charToRaw(paste(expr,collapse='\n'))
  l <- length(r) + 1L
  writeBin(as.integer(c(0x42L, l + 4L,0L ,0L ,4L + l * 256L)), c, endian="little")
  writeBin(r, c)
  writeBin(raw(1), c)
  b <- readBin(c, "int", 4, endian="little")
  if (!length(b)) { close(c); stop("Rserve connection timed out and closed") }
  msgLen <- b[2]
  if (msgLen > 0) a <- readBin(c,"raw",msgLen)
  if (b[1]%%256 != 1) stop("RSserverEval failed with error: ",b[1]%/%0x1000000)
  invisible(b[1]%%256 == 1)  
}

RSserverSource <- function(c, file) {
  if (!is.character(file) || length(file) != 1) stop("`file' must be a string")
  r <- charToRaw(file)
  l <- length(r) + 1L
  writeBin(as.integer(c(0x45L, l + 4L,0L ,0L ,4L + l * 256L)), c, endian="little")
  writeBin(r, c)
  writeBin(raw(1), c)
  b <- readBin(c, "int", 4, endian="little")
  if (!length(b)) { close(c); stop("Rserve connection timed out and closed") }
  msgLen <- b[2]
  if (msgLen > 0) a <- readBin(c,"raw",msgLen)
  if (b[1]%%256 != 1) stop("RSserverSource failed with error: ",b[1]%/%0x1000000)
  invisible(b[1]%%256 == 1)  
}

RSshutdown <- function(c, pwd=NULL, ctrl=FALSE) {
  if (ctrl) {
    writeBin(c(0x44L, 0L, 0L, 0L), c, endian="little")
    b <- readBin(c, "int", 4, endian="little")
    if (!length(b)) { close(c); stop("Rserve connection timed out and closed") }
    msgLen <- b[2]
    if (msgLen > 0) a <- readBin(c,"raw",msgLen)
    if (b[1]%%256 != 1) stop("ctrlShutdown failed with error: ",b[1]%/%0x1000000)
    invisible(b[1]%%256 == 1)  
  } else {
    # FIXME: we ignore pwd and don't check error status
    writeBin(as.integer(c(4, 0, 0, 0)), c, endian="little")
  }
}
