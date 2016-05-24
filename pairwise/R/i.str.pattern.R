str.pattern<-function (x,d=1)
{
  # converts x into string patterns by row == default, or by colum (d = 2)
  # NAs in the (data)matrix x result in additional response patterns
  # func. by joerg-henrik heine jhheine(at)googlemail.com
ppaste <- function(z){paste(z, collapse = "")} # collaps funktion pattern strings
strpat<-(apply(x, d, ppaste)) # creates pattern strings
return(strpat)
}