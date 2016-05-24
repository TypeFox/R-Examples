#################################################################
# 
# File:         read.cdb.R
# Purpose:      Read a cdb database
#
# Created:      20130416
# Authors:      Emilio Torres Manzanera
#
# Modifications: 
#
#################################################################

read.cdb <- function(file, type=c("cdb","txt")) {
  data <- data.frame(key=as.character(NULL), value=as.character(NULL))
  if(type[1] == "cdb" ) {
    cdbunpack <- function( buffer ) {
      if( class(buffer) != "raw")
        stop("Not raw data")
      if( length(buffer) != 4)
        stop("No length 4")
      require("bitops")
      n <- buffer[4]
      n <- bitShiftL(n, 8)
      n <- bitOr(n, buffer[3])
      n <- bitShiftL(n, 8)
      n <- bitOr(n, buffer[2])
      n <- bitShiftL(n, 8)
      n <- bitOr(n, buffer[1])
      n
    }
    zz <- file(file,"rb")
    pointerend <- cdbunpack(readBin(zz, "raw", n=4, size= 1, signed = FALSE))
    dummy <- readChar(zz, 2048 - 4, useBytes = TRUE)
    pos <- 2048
    data <- data.frame(key=as.character(NULL), value=as.character(NULL))
    while( pos < pointerend ) {
      klen <- cdbunpack(readBin(zz, "raw", n=4, size= 1, signed = FALSE))
      vlen <- cdbunpack(readBin(zz, "raw", n=4, size= 1, signed = FALSE))
      key <-  readChar(zz,klen, useBytes = TRUE)
      value <- readChar(zz,vlen, useBytes = TRUE)
      pos <- pos + 8 + klen + vlen
      data <- rbind(data, cbind(key,value))
    }
    close(zz)
  } else {
    readalineofcdbtxtformat <- function(zz) {
      dummy <-  readChar(zz,1) # Read the "+"
      if(length(dummy)!=1 | dummy != "+" )
        return(NULL)
      klenstring <- ""
      while( (z <- readChar(zz,1,useBytes=TRUE)) != ",") {
        klenstring <- paste(klenstring,z,sep="")}
      vlenstring <- ""
      while( (z <- readChar(zz,1,useBytes=TRUE)) != ":") {
        vlenstring <- paste(vlenstring,z,sep="")}
      key <- readChar(zz,strtoi(klenstring),useBytes=TRUE)
      dummy <- readChar(zz,2) # Read the "->"
      value <- readChar(zz,strtoi(vlenstring),useBytes=TRUE)
      dummy <- readChar(zz,1) # Read the "\n"
      if( strtoi(klenstring) != nchar(key, type = "bytes") |
         strtoi(vlenstring) != nchar(value, type = "bytes") )
        stop("Error in length of key or value!")
      cbind(key,value)
    }
    zz <- file(file,"r")
    if(!isOpen(zz))
      stop(paste("File ", file, " does not exist"))
    while( length(register <- readalineofcdbtxtformat(zz)) == 2 ) {
      data <- rbind(data, register)
    }
    close(zz)
  }
  data
}
