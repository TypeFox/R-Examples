words <- function(length = 3, alphabet = s2c("acgt") )
{
  if( length == 1 )
    return( alphabet )
  else
    kronecker( alphabet, words(length - 1, alphabet ), paste, sep = "")
}
