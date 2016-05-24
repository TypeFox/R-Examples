vecprint <- function(x, n.char=78) {

################################################################################
# Function: vecprint
# Programmer: Tom Kincaid
# Date: February 6, 2004
# Last Revised:  May 3, 2004
# Description:
#   This function takes an input vector and outputs a character string with line
#   breaks inserted so that, whenever possible, no line in the string exceeds
#   the input value n.char, which is set to 78 characters by default.  The input
#   vector is coereced to mode character.  When an element of the input vector
#   is greater than n.char characters in length, then that element is inserted
#   in the output character string as an individual line. 
#   Input:
#      x = a vector.
#      n.char = the maximum number of characters per line.  The default is 78.
#   Output:
#      A character string that is suitable for printing by the functions: stop,
#      warning, or cat.
#   Other Functions Required: None
#   Examples:
#      sites <- paste("Site Number", 1:50)
#      sites.str <- vecprint(sites)
#      cat(sites.str)
#
#      temp <- c(1, 5, 21:25, 33:37)
#      sites.str <- vecprint(sites[temp])
#      warning(paste("\nThe following site ID values were removed from the
#         analysis:\n", sites.str, sep=""))
################################################################################

   x <- as.character(x)
   n <- length(x)
   nc <- nchar(x)

   i <- 1
   j <- 1
   nc.sum <- 0
   x.str <- character()
   while(j <= n) {
      if((nc.sum + nc[j]) <= n.char) {
         nc.sum <- nc.sum + nc[j] + 2
         j <- j+1
      } else {
         if(i == j) {
            x.str <- paste(x.str, x[j], "\n")
         } else {
            j <- j-1
            x.str <- paste(x.str, paste(x[i:j], collapse=", "), "\n")
         }
         j <- j+1
         i <- j
         nc.sum <- 0
      }
   }
   if(i < j) {
      j <- j-1
      x.str <- paste(x.str, paste(x[i:j], collapse=", "), "\n")
   }

   x.str
}

