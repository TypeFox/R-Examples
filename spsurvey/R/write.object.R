write.object <- function(obj, n.digits=2, r.names=TRUE, c.names=TRUE, r.cex=1,
   c.cex=1, miss="NA") {

################################################################################
# Function: write.object
# Purpose: Write contents of an object to a plot
# Programmer: Tom Kincaid
# Date: February 5, 2001
# Last Revised: May 7, 2010
# Description:
#   This function writes the contents of an object to a plot.  The object may be
#   either a data frame or a matrix.  Values in the input data frame or matrix
#   must be of class numeric, character, or factor.
# Arguments:
#   obj = the object (either a data frame or a matrix).
#   n.digits = number of digits after the decimal point for numeric values.  The
#     default is 2.
#   r.names = a logical value that indicates whether to print the row names,
#     where TRUE = print the row names and  FALSE = do not print the row names.
#     The default is TRUE.
#   c.names = a logical value that indicates whether to print the column names,
#     where TRUE = print the column names and  FALSE = do not print the column
#     names.  The default is TRUE.
#   r.cex = character expansion parameter for the row labels.  The default is 1.
#   c.cex = character expansion parameter for the column labels.  The default is
#     1.
#   miss = the missing value code expressed as a character string.  The default
#    is "NA".
# Results:
#   The function returns NULL.  Side effect of the function is to write contents
#   of the input object to a plot.
# Other Functions Required:
#   input.format - format an input value
# Examples:
#   z <- rnorm(100)
#   z.mean <- c(tapply(z, rep(1:4, rep(25,4)), mean), mean(z))
#   z.sd <- sqrt(c(tapply(z, rep(1:4, rep(25,4)), var), var(z)))
#   z.upper <- z.mean+1.96*z.sd
#   z.lower <- z.mean-1.96*z.sd
#   obj <- data.frame(rbind(z.mean, z.sd, z.upper, z.lower))
#   dimnames(obj) <- list(c("Mean Estimate", "Standard Deviation", "Lower 95\%
#     Conf. Bound", "Upper 95\% Conf. Bound"), c(paste("Stratum", 1:4, sep=""),
#     "AllStrata"))
#   write.object(obj, n.digits=3, r.cex=0.75)
#
#   obj <- data.frame(matrix(round(5 + runif(30), 1), nrow=6))
#   colnames(obj) <- c("United States", "Russia", "Germany", "Japan", "France")
#   write.object(obj, n.digits=1, r.names=FALSE)
################################################################################

# If the object is a matrix, convert to a data frame

   if(!is.data.frame(obj)) obj <- data.frame(obj)

# Assign the number of rows and number of columns for the plot

   dim.obj <- dim(obj)
   if (r.names)
      n.col <- dim.obj[2] + 1
   else
      n.col <- dim.obj[2]
   if (c.names)
      n.row <- dim.obj[1] + 1
   else
      n.row <- dim.obj[1]

# Create the plot area

   plot(seq(n.col), seq(n.col), xlim=c(0.5, n.col+0.5), ylim=c(0.5, n.row+0.5), type="n", axes=FALSE, xlab="", ylab="")

# Plot the object

   if (r.names) {
      if (c.names) {
         text(0.5, rev(1:(n.row-1)), dimnames(obj)[[1]], adj=0, cex=r.cex)
         text(2:n.col, n.row, dimnames(obj)[[2]], adj=0.5, cex=c.cex)
         for (i in 1:(n.row-1)) {
            for(j in 1:(n.col-1)) {
               text(j+1, n.row-i, input.format(obj[i,j], n.digits, miss), adj=0.5)
            }
         }
      } else {
         text(0.5, rev(1:n.row), dimnames(obj)[[1]], adj=0, cex=r.cex)
         for (i in 1:n.row) {
            for(j in 1:(n.col-1)) {
               text(j+1, n.row-i+1, input.format(obj[i,j], n.digits, miss), adj=0.5)
            }
         }
      }
   } else {
      if (c.names) {
         text(seq(n.col), n.row, dimnames(obj)[[2]], adj=0.5, cex=c.cex)
         for (i in 1:(n.row-1)) {
            for(j in 1:n.col) {
               text(j, n.row-i, input.format(obj[i,j], n.digits, miss), adj=0.5)
            }
         }
      } else {
         for (i in 1:n.row) {
            for(j in 1:n.col) {
               text(j, n.row-i+1, input.format(obj[i,j], n.digits, miss), adj=0.5)
            }
         }
      }
   }

# Place a box around the plot area

   box("plot", lwd=2)

# Return a NULL object

   invisible(NULL)
}

