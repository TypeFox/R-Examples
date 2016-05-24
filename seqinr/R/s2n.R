##############################################################################
#
# simple numerical encoding of a DNA sequence that by default
# is independent of locale.
#
##############################################################################

s2n <- function(seq, levels = s2c("acgt"), base4 = TRUE, forceToLower = TRUE)
{
  #
  # Check that sequence is a vector of chars:
  #
  if(nchar(seq[1]) > 1) stop("sequence is not a vector of chars")
  #
  # Force to lower-case letters if requested:
  #
  if(forceToLower) seq <- tolower(seq)

  if( base4 )
    unclass(factor(seq, levels = levels ) ) - 1
  else
    unclass(factor(seq, levels = levels ) )
}