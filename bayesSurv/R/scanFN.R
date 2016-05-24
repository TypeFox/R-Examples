#*** scanFN.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##           CREATED:  01/02/2007
##
## PURPOSE: Customized 'scan'
##
## FUNCTIONS: scanFN
## 
#* ********************************************************************************* */
scanFN <- function(file, quiet=FALSE)
{
  xname <- scan(file, what=character(0), nlines=1, quiet=quiet)
  nx <- length(xname)
  x <- matrix(scan(file, what=double(0), skip=1, quiet=quiet), ncol=nx, byrow=TRUE)
  colnames(x) <- xname
  return(as.data.frame(x))  
}

