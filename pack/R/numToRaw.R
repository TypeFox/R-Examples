#-------------------------------------------------------------------------#
# pack R package, copyright (C) Joshua M. Ulrich, 2007-2008               #
# Distributed under GNU GPL version 3                                     #
#-------------------------------------------------------------------------#

'numToRaw' <-
function(x, nBytes=1) {

  # from 'wle' package
  # Supporting function to convert numbers to a
  # raw vector of length 'nBytes'

  pos <- ifelse( x==0, 1, floor(log(x,2))+1 )
  if( pos <= nBytes*8 )
    pos <- nBytes*8
  else
    stop('the value of \"nBytes\" is too small')

  bin <- rep(0, pos)
  for( i in pos:1 ) {
    bin[i] <- floor(x/2^(i-1))
    x <- x-((2^(i-1))*bin[i])
  }
  return(packBits(as.integer(bin)))
}
