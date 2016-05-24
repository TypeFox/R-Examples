###########################################################################
#  Internal functions (not documented)
############################################################################

"valid.one" <- function(x,type)
{
  if(!eval(call(paste("is.",type,sep=""),x))||length(x)!=1) {
    error <- paste("'",deparse(substitute(x)),"' must be a length-one object of type ",type,sep="")
    stop(error)
  } 
}

