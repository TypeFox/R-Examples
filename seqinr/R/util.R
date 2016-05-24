########################
# char to string
########################

c2s <- function( chars = c("m","e","r","g","e","d") )
{
  return( paste( chars, collapse = "" ) )
}

###########################
# string to char
############################

s2c <- function (string) 
{
  if(is.character(string) && length(string) == 1){
    return(.Call("s2c", string, PACKAGE = "seqinr"))
  } else {
    warning("Wrong argument type in s2c(), NA returned")
    return(NA)
  }
}



###########################
# Conversion of the numeric encoding of a DNA sequence into
# a vector of chars
############################

n2s <- function(nseq, levels = c("a", "c", "g", "t"), base4 = TRUE)
{
  if( base4 )
    levels[nseq + 1]
  else
    levels[nseq]
}



##########################################
# Convert one-letter code to 3-letters code for amino-acids
##########################################

aaa <- function( aa )
{
  aa3 <- c("Stp", "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile",
           "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr",
           "Val", "Trp", "Tyr") # One letter code order
  if(missing(aa)) return(aa3)
  aa1 <- a()

  convert <- function( x )
  {
    if( all( x != aa1 ) )
    { 
      warning("Unknown one letter code for aminoacid")
      return( NA )
    }
    else
    {
      return( aa3[which( x == aa1 )] )
    }
  }
  return( as.vector(unlist(sapply( aa, convert ) ) ) )
}

##########################################
# Conversion 3-letters code to one letter code for amino-acids
##########################################

a <- function( aa )
{
  aa1 <- s2c("*ACDEFGHIKLMNPQRSTVWY")
  if(missing(aa)) return(aa1)
  aa3 <- aaa()

  convert <- function( x )
  {
    if( all( x != aa3 ) )
    { 
      warning("Unknown 3-letters code for aminoacid")
      return( NA )
    }
    else
    {
      return( aa1[which( x == aa3 )] )
    }
  }
  return( as.vector(unlist(sapply( aa, convert ) ) ) )
}


