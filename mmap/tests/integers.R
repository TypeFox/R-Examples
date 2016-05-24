library(mmap)
tmp <- tempfile()

##### int8() #####
# write binary from min to max signed 2^8
test.int8 <- function() {
  writeBin(-128:127L, tmp, size=1)
  
  m <- mmap(tmp, int8())  # signed 1 byte integers
  if( !all(m[] == (-128:127L)) )
    stop("m[] == (-128:127L)")
  
  # test replacement
  m[] <- -128L
  if( !all(m[] == -128))
    stop("m[] == -128")
  munmap(m)
}




#### uint8() ####
test.uint8 <- function(on=TRUE) {
  if( !isTRUE(on)) {
    cat("test.uint8 disabled\n")
    return(NULL)
  }
  writeBin(0:255L, tmp, size=1)
  m <- mmap(tmp, uint8())  # unsigned 1 byte integers
  if( !all(m[] == 0:255L) )
    stop("m[] == 0:255L")
  
  # test replacement
  m[] <- 255L;
  if( !all(m[] == 255L))
    stop("m[] == 255L")
  munmap(m)
}




#### int16() ####
test.int16 <- function(on=TRUE) {
  if( !isTRUE(on)) {
    cat("test.int16 disabled\n")
    return(NULL)
  }
  writeBin(-32768:32767L, tmp, size=2)
  m <- mmap(tmp, int16())  # signed 2 byte integers
  if( !all(m[] == -32768:32767L) )
    stop("m[] == -32768:32767L")
  
  # test replacement
  m[] <- -32768L
  if( !all(m[] == -32768L))
    stop("m[] == -32768L")
  munmap(m)
}




#### uint16() ####
test.uint16 <- function(on=TRUE) {
  cat("checking test.uint16...")
  if( !isTRUE(on)) {
    cat("test.uint16 disabled\n")
    return(NULL)
  }
  writeBin(0:65535L, tmp, size=2)
  m <- mmap(tmp, uint16())  # unsigned 2 byte integers
  if( !all(m[] == 0:65535L) )
    stop("m[] == 0:65535L")
  
  # test replacement
  m[] <- 65535L
  if( !all(m[] == 65535L))
    stop("m[] == 65535L")
  munmap(m)
  cat("OK\n")
}




#### int24() ####
test.int24 <- function(on=TRUE) {
  cat("checking test.int24...")
  if( !isTRUE(on)) {
    cat("test.int24 disabled\n")
    return(NULL)
  }
  ints <- as.integer(seq(-8388607L,8388607L,length.out=11))
  writeBin(rep(as.raw(0),33), tmp)
  m <- mmap(tmp, int24())  # signed 3 byte integers
  m[] <- ints
  if( !all(m[] == ints) )
    stop("m[] == ints")
  munmap(m)
  cat("OK\n")
}




#### uint24() ####
test.uint24 <- function(on=TRUE) {
  cat("checking test.uint24...")
  if( !isTRUE(on)) {
    cat("test.uint24 disabled\n")
    return(NULL)
  }
  ints <- as.integer(seq(0,16777215L,length.out=100))
  writeBin(rep(as.raw(0),300), tmp)
  m <- mmap(tmp, uint24())  # signed 3 byte integers
  m[] <- ints
  if( !all(m[] == ints) )
    stop("m[] == ints")
  munmap(m)
  cat("OK\n")
}




#### int32() ####
test.int32 <- function(on=TRUE) {
  cat("checking test.int32...")
  if( !isTRUE(on)) {
    cat("test.int32 disabled\n")
    return(NULL)
  }
  writeBin(-1e6:1e6L, tmp, size=4)
  m <- mmap(tmp, int32())  # unsigned 2 byte integers
  if( !all(m[] == -1e6:1e6L) )
    stop("m[] == -1e6:1e6L")
  
  # test replacement
  m[] <- .Machine$integer.max
  if( !all(m[] == .Machine$integer.max))
    stop("m[] == .Machine$integer.max")
  munmap(m)
  cat("OK\n")
}




#### int64() ####
test.int64 <- function(on=TRUE) {
  cat("checking test.int64...")
  if( !isTRUE(on)) {
    cat("test.int32 disabled\n")
    return(NULL)
  }
  writeBin(0.0, tmp)
  m <- mmap(tmp, int64())  # signed 8 byte integers as doubles
  m[] <- 2^40
  if( !all(m[] == 2^40) )
    stop("m[] == 2^40")
  munmap(m)
  cat("OK\n")
}

test.int8()
test.uint8()
test.int16()
test.uint16()
test.int24()
test.uint24()
test.int32()
test.int64(FALSE)
