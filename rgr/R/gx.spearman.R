gx.spearman <-
function(xx, ifclr = FALSE, ifwarn = TRUE)
{
     # Function to compute Spearman correlation coefficients and their
     # significance for a matrix, rejecting any rows containing NAs,
     # and display the results to 3 significant figures.
     #
     # NOTE: Prior to using this function the data frame/matrix containing the
     # variables, xx, must be run through ltdl.fix.df to convert any <dl -ve
     # values to positive half that value, and set zero2na = TRUE if it is
     # required to convert any zero values or other numeric codes representing 
     # blanks to NAs.
     #
     if(!is.matrix(xx)) stop(deparse(substitute(xx)), " is not a Matrix")
     # Remove any vectors containing NAs
     temp.x <- remove.na(xx)
     x <- temp.x$x
     if(ifclr) {
         x <- clr(x, ifwarn = ifwarn)
         cat("Data have been Centred Log-Ratio transformed\n")
     }
     #Compute correlation matrix
     r <- cor(x, method = "spearman")
     df.t <- temp.x$n - 2
     df.term <- sqrt(df.t)
     for(i in 2:temp.x$m) {
         for(j in 1:(i - 1))
             r[i, j] <- pt((abs(r[i, j]) * df.term)/
                        sqrt(1 - r[i, j] * r[i, j]), df.t)
     }
     r <- round(r, 3)
     for(i in 1:temp.x$m) r[i,i] <- NA
     cat("Spearman Correlation Coefficients and their Statistical Significance,",
         "\nupper and lower triangles, respectively,",
         paste("for matrix ", deparse(substitute(xx)), ", N = ", temp.x$n, "\n\n", sep = ""))
     print(r, na.print = " ")
     cat("\n")
     invisible()
}
