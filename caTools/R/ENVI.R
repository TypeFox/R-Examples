#===========================================================================#
# caTools - R library                                                       #
# Copyright (C) 2005 Jarek Tuszynski                                        #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

write.ENVI = function(X, filename, interleave=c("bsq", "bil", "bip") ) 
{ # write matrix or data cube to binary ENVI file
  if (is.vector(X)) {
    nCol = length(X)
    nRow <- nBand <- 1
  } else {
    d = dim(X)
    nRow  = d[1]
    nCol  = d[2]
    nBand = prod(d)/(nRow*nCol)
  }
  dim(X) = c(nRow, nCol, nBand)    # make it into 3D array in case it was not
  
  # check data type
  data.type = 0
  if (is.double (X)) data.type = 5 # 64-bit double
  if (is.integer(X)) data.type = 3 # 32-bit int
  if (is.complex(X)) data.type = 9 # 2x64-bit complex<double>
  if (data.type == 0) {            # do not know what is it -> make it a double
    X = as.double(X) 
    data.type = 5 
  } 
  
  # change interleave and store tha data
  interleave = match.arg(interleave)
  if      (interleave=="bil") X=aperm(X, c(2,3,1))  # R's [row,col,band] -> bil [col,band,row] 
  else if (interleave=="bip") X=aperm(X, c(3,2,1))  # R's [row,col,band] -> bip [band,col,row] 
  else if (interleave=="bsq") X=aperm(X, c(2,1,3))  # R's [row,col,band] -> bsq [col,row,band] 
  writeBin(as.vector(X), filename)                  # write Envi file

  # write header file
  out  = "ENVI\ndescription = { R-language data }\n"
  out  = paste(out, "samples = ", nCol, "\n", sep="")
  out  = paste(out, "lines = ", nRow, "\n", sep="")
  out  = paste(out, "bands = ", nBand, "\n", sep="")
  out  = paste(out, "data type = ",data.type,"\n", sep="")
  out  = paste(out, "header offset = 0\n", sep="")
  out  = paste(out, "interleave = ",interleave,"\n", sep="")   # interleave is assumed to be bsq - in case of 1 band images all 3 formats are the same 
  ieee = if(.Platform$endian=="big") 1 else 0       # does this machine uses ieee (UNIX) format? or is it intel format?
  out  = paste(out, "byte order = ", ieee, "\n", sep="")
  cat(out, file=paste(filename, ".hdr", sep=""))
  invisible(NULL)
}

# =======================================================================================

read.ENVI = function(filename, headerfile=paste(filename, ".hdr", sep=""))  
{  # read matrix or data cube from binary ENVI file
  
  # parse header file 
  nCol <- nRow <- nBand <- data.type <- header.offset <- byte.order <- (-1)
  interleave = "bsq"
  if (!file.exists(headerfile)) stop("read.ENVI: Could not open input header file: ", headerfile)
  Lines  = read.table(headerfile, sep="=", strip.white=TRUE, row.names = NULL, as.is=TRUE, fill=TRUE)
  Fields = c("samples", "lines", "bands", "data type", "header offset", "interleave", "byte order")
  for (i in 1:nrow(Lines)) {
    Lab = tolower(Lines[i,1])
    Lab = gsub("[ ]+", " ", Lab) # Replace all multiple spaces with a single space 
    j = match(Lab, Fields)
    Val = Lines[i,2]
    if (length(j) == 1)
     switch( j, 
       nCol          <- as.integer(Val),
       nRow          <- as.integer(Val),
       nBand         <- as.integer(Val),
       data.type     <- as.integer(Val),
       header.offset <- as.integer(Val),
       interleave    <- gsub(" ", "", Val),
       byte.order    <- as.integer(Val)
     )
   }

  if (nCol <= 0 | nRow <= 0 | nBand <= 0) 
    stop("read.ENVI: Error in input header file ", headerfile, " data sizes missing or incorrect", nRow, nCol, nBand)
  if (! ( data.type %in% c(1,2,3,4,5,9,12) ) ) 
    stop("read.ENVI: Error in input header file ", headerfile, " data type is missing, incorrect or unsupported ")
  
  # read the data binary file
  ieee = if(.Platform$endian=="big") 1 else 0       # does this machine uses ieee (UNIX) format? or is it intel format?
  endian = if(ieee==byte.order | byte.order<0) .Platform$endian else "swap"  
  size   = nRow*nCol*nBand
  if (!file.exists(filename)) stop("read.ENVI: Could not open input file: ", filename)
  f = file(filename, "rb")
  if (header.offset>0) readBin(f, raw(), n=header.offset)
  switch( data.type, 
    X <- readBin(f, integer(), n=size, size=1, signed=FALSE),  # data.type==1 -> 1-byte unsigned integer (char)
    X <- readBin(f, integer(), n=size, size=2, endian=endian), # data.type==2 -> 2-byte short 
    X <- readBin(f, integer(), n=size, endian=endian),         # data.type==3 -> 4-byte int
    X <- readBin(f, double() , n=size, size=4, endian=endian), # data.type==4 -> 4-byte float 
    X <- readBin(f, double() , n=size, endian=endian), , , ,   # data.type==5 -> 8-byte double
    X <- readBin(f, complex(), n=size, endian=endian), , ,     # data.type==9 -> 2x8-byte complex<double>
    X <- readBin(f, integer(), n=size, size=2, endian=endian, signed=FALSE) # data.type==12 -> 2-byte unsigned short integer
  )
  close(f)

  Fields = c("bil", "bip", "bsq")
  j = match(interleave, Fields)
  if (length(j)==0) stop("read.ENVI: Error in input header file ", headerfile, " incorrect interleave type")
  switch(j,        
   { dim(X)<-c(nCol,nBand,nRow); X<-aperm(X, c(3,1,2)); }, # bil [col,band,row] -> R's [row,col,band]
   { dim(X)<-c(nBand,nCol,nRow); X<-aperm(X, c(3,2,1)); }, # bip [band,col,row] -> R's [row,col,band]
   { dim(X)<-c(nCol,nRow,nBand); X<-aperm(X, c(2,1,3)); }  # bsq [col,row,band] -> R's [row,col,band]
  )
  if (nBand==1) dim(X)=c(nRow, nCol)
  return(X)
} 

