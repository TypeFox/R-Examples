### <======================================================================>
#
# This is a helper function for ecd constructor
# It extends is.numeric for mpfr
#
"is.numericMpfr" <- function(x)
{
    ifelse( is.numeric(x) | class(x)=="mpfr" | class(x)=="numericMpfr", TRUE, FALSE ) 
}
### <---------------------------------------------------------------------->
