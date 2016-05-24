#-------------------------------------------------------------------------#
# pack R package, copyright (C) Joshua M. Ulrich, 2007-2008               #
# Distributed under GNU GPL version 3                                     #
#-------------------------------------------------------------------------#

'rawToNum' <-
function(x, nBytes=1) {
  
  # Supporting function to convert raw vectors
  # to numbers according to 'nBytes'

  i <- as.logical(rawToBits(x))
  num <- sum(2^.subset(0:(nBytes*8-1), i))

  return(num)
}
