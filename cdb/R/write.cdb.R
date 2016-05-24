#################################################################
# 
# File:         write.cdb.R
# Purpose:      write a data frame in a text file with a cdb record format
#
# Created:      20130416
# Authors:      etm
#
# Modifications: 
#
#################################################################

write.cdb <-  function(x, file, type=c("cdb","txt")) {
  if (!is.data.frame(x)) 
    x <- data.frame(x)
  if( dim(x)[2]!= 2 )
    stop("x must be pairs of (key,value) registers")
  names(x) <- c("key","val")
  if(type[1] == "cdb" ) {
    message("Not yet implemented. Use the txt option.")
  }
  else {
  x$key <- as.character(x$key)
  x$val <- as.character(x$val)
  x$klen <- nchar(x$key, type = "bytes")
  x$vlen <- nchar(x$val, type = "bytes")
  text <- paste("+",x$klen,",",x$vlen,":",x$key,"->",x$val,sep="")
  fileConn<-file(file)
  ## We have to add a empty line at the end of the file
  writeLines(c(text,"\n"), fileConn,sep="\n")
  close(fileConn)
  }
}



