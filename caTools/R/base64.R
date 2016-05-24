#===========================================================================#
# caTools - R library                                                       #
# Copyright (C) 2005 Jarek Tuszynski                                        #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

#===============================================================================
# The Base64 encoding is designed to encode arbitrary binary information for 
# transmission by electronic mail. It is defined by MIME (Multipurpose Internet 
# Mail Extensions) specification RFC 1341, RFC 1421, RFC 2045 and others. 
# Triplets of 8-bit octets are encoded as groups of four characters, each 
# representing 6 bits of the source 24 bits. Only a 65-character subset 
# ([A-Z,a-z,0-9,+,/,=]) present in all variants of ASCII and EBCDIC is used, 
# enabling 6 bits to be represented per printable character
#===============================================================================

base64encode = function(x, size=NA, endian=.Platform$endian)
{
   if ((typeof(x)!="character")&(typeof(x)!="raw")) x = writeBin(x, raw(), size=size, endian=endian)
   if ((typeof(x)=="character")&(typeof(x)!="raw")) {nlen<- nchar(x);x = writeBin(x, raw(), size=size, endian=endian);length(x)<- nlen}
   x = as.integer(x)
   ndByte = length(x)            # number of decoded bytes
   nBlock = ceiling(ndByte / 3)  # number of blocks/groups
   neByte = 4 * nBlock           # number of encoded bytes

   # add padding if necessary, to make the length of x a multiple of 3
   if (ndByte < 3*nBlock) x[(ndByte+1) : (3*nBlock)] = 0;
   dim(x) = c(3, nBlock)         # reshape the data
   y = matrix(as.integer(0), 4, nBlock)  # for the encoded data
   
   #-------------------------------------------
   # Split up every 3 bytes into 4 pieces
   #   x = aaaaaabb bbbbcccc ccdddddd
   # to form
   #   y = 00aaaaaa 00bbbbbb 00cccccc 00dddddd
   # than convert y to integers in 0-63 range
   # This section is based on Matlab code by Peter Acklam
   # http://home.online.no/~pjacklam/matlab/software/util/datautil/
   #-------------------------------------------
   y[1,] = bitShiftR(x[1,], 2) # 6 highest bits of x(1,:)
   y[2,] = bitOr(bitShiftL(x[1,], 4), bitShiftR(x[2,], 4))
   y[3,] = bitOr(bitShiftL(x[2,], 2), bitShiftR(x[3,], 6))
   y[4,] = x[3,]
   y = bitAnd(y, 63)           # trim numbers to lower 6-bits
   
   #----------------------------------
   # Perform the following mapping
   #   0  - 25  ->  A-Z
   #   26 - 51  ->  a-z
   #   52 - 61  ->  0-9
   #   62       ->  +
   #   63       ->  /
   #----------------------------------
   alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
   alpha = strsplit(alpha, NULL)[[1]]  # convert string to array of characters
   z = alpha[y+1]                      # rearrange characters 

   #-------------------------
   # Add padding if necessary.
   #-------------------------
   npbytes = 3 * nBlock - ndByte   # number of padding bytes needed
   if (npbytes>0) z[(neByte-npbytes+1) : neByte] = '='  # '=' is used for padding
   z = paste(z, collapse = "")       # combine characters into a string
   return (z)
}

#====================================================================

base64decode = function(z, what, size=NA, signed = TRUE, endian=.Platform$endian)
{  
  if (!is.character(z)) 
    stop("base64decode: Input argument 'z' is suppose to be a string")
  if (length(z)==1) z = strsplit(z, NULL)[[1]] # convert string to array of characters
  if (length(z)%%4!=0) 
   warning("In base64decode: Length of base64 data (z) not a multiple of 4.")
  #-----------------------------------
  # Now perform the following mapping
  #   A-Z  ->  0  - 25
  #   a-z  ->  26 - 51
  #   0-9  ->  52 - 61
  #   +    ->  62
  #   /    ->  63
  #   =    ->  64  - special padding character
  #  otherwise -1
  #-----------------------------------
  alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="
  alpha = strsplit(alpha, NULL)[[1]]    # convert string to array of characters
  y     = match(z, alpha, nomatch=-1)-1 # lookup number of each character
  if (any(y == -1)) 
    stop("base64decode: Input string is not in Base64 format")
  if (any(y == 64)) y = y[y != 64]      # remove padding
  neByte = length(y);                   # number of encoded bytes
  nBlock = ceiling(neByte/4);           # number of blocks/groups
  ndByte = 3 * nBlock                   # number of decoded bytes
  
  # add padding if necessary
  if (neByte < 4*nBlock) y[(neByte+1) : (4*nBlock)] = 0;
  dim(y) = c(4, nBlock);                # shape into a matrix
  x = matrix(as.integer(0), 3, nBlock); # for the decoded data
 
  #---------------------------------------------
  # Rearrange every 4 bytes into 3 bytes
  #    y = 00aaaaaa 00bbbbbb 00cccccc 00dddddd
  # to form
  #    x = aaaaaabb bbbbcccc ccdddddd
  # This section is based on Matlab code by Peter Acklam
  # http://home.online.no/~pjacklam/matlab/software/util/datautil/
  #---------------------------------------------
  x[1,] = bitOr(bitShiftL(y[1,], 2), bitShiftR(y[2,], 4))
  x[2,] = bitOr(bitShiftL(y[2,], 4), bitShiftR(y[3,], 2))
  x[3,] = bitOr(bitShiftL(y[3,], 6), y[4,])
  x = bitAnd(x, 255) # trim numbers to lower 8-bits
  
  # remove padding
  if (neByte %% 4 == 2) x = x[1:(ndByte-2)]
  if (neByte %% 4 == 3) x = x[1:(ndByte-1)]
  
  # perform final conversion from 'raw' to type given by 'what'
  r = as.raw(x)
  TypeList = c("logical", "integer", "double", "complex", "character", "raw", 
               "numeric", "int")
  if (!is.character(what) || length(what) != 1 || !(what %in% TypeList)) 
    what <- typeof(what)
  if (what=="raw") return(r)
  if (is.na(size)) size = switch(match(what, TypeList), 4, 4, 8, 16, 2, 1, 8, 4)
  if (what=="character") {rlen<- size*ceiling(length(r)/size);length(r)<- rlen}  
  n = length(r)
  if (n%%size) stop("raw2bin: number of elements in 'r' is not multiple of 'size'")
  x = readBin(r, what, n = n%/%size, size=size, signed=signed, endian=endian)
  if (what=="character")  x = paste(x, collapse = "") # convert arrays of characters to strings
  return (x)
}

