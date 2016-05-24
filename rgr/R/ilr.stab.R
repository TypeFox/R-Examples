ilr.stab <-
function(xx1, xx2, ifwarn = T) 
{
     # Function to compute the robust ilr stability for 2 parts of a composition
     # using the squared MAD in place of the variance.
     #
     # NOTE: Prior to using this function the data frame/matrix containing the
     # variables, xx1 and xx2, must be run through ltdl.fix.df to convert any <dl
     # -ve values to positive half that value, and set zero2na = TRUE if it is
     # required to convert any zero values or other numeric codes representing 
     # blanks to NAs.
     #
     if (length(xx1) != length(xx2)) 
         stop("Input vectors must be of equal length\n")
     if(ifwarn) cat("  ** Are the data in the same measurement units? **\n" )
     temp.x <- remove.na(cbind(xx1, xx2))
     x1 <- temp.x$x[1:temp.x$n, 1]
     x2 <- temp.x$x[1:temp.x$n, 2]
     ilr.mad <- mad(log(x1 / x2)/1.4142)
     ilr.stab <- exp(-ilr.mad * ilr.mad)
     return(ilr.stab = ilr.stab)
}
