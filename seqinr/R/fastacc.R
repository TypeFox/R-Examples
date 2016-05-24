#                 fastacc: Fast Allele in Common Count
#
# Computes the number of common alleles between a target and a database
#
fastacc <- function(target, database)
{
  #
  # Check arguments:
  #
  if(!is.raw(target)) stop("raw vector expected for target")
  noc <- length(target)
  if(noc < 1) stop("empty target")
  if(!is.raw(database)) stop("raw vector expected for database")
  if(length(database) %% noc != 0) stop("database length must be a multiple of target length")
  #
  # n is the total number of entries in DB
  #
  n <- length(database)/noc
  n <- as.integer(n)
  #
  # Pre-compute the table giving the number of bits per octet
  #
  bits_in_char <- sapply(as.raw(0:255), function(x) sum(as.integer(rawToBits(x))))
  bits_in_char <- as.integer(bits_in_char) # just to make sure
  
  res <- .Call("fastacc", bits_in_char, target, database, noc, n, PACKAGE = "seqinr")
  
  return(res)
}
