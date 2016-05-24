gx.vm <-
function(xx, ifwarn = TRUE)
{
     # Function to compute an Aitchison variability matrix, rejecting any of the
     # input matrix containing NAs, and displaying the results to 3 significant
     # figures.
     #
     # NOTE: Prior to using this function the data frame/matrix containing the
     # variables, xx, must be run through ltdl.fix.df to convert any <dl -ve
     # values to positive half that value, and set zero2na = TRUE if it is
     # required to convert any zero values or other numeric codes representing 
     # blanks to NAs.
     #
     if(!is.matrix(xx)) stop(deparse(substitute(xx)), " is not a Matrix")
     if(ifwarn) cat("  ** Are the data all in the same measurement units? **\n" )
     # Remove any vectors containing NAs
     temp.x <- remove.na(xx)
     x <- temp.x$x
     nx <- temp.x$m
     r <- matrix(nrow = nx, ncol = nx)
     dimnames(r)[[2]] <- dimnames(r)[[1]] <- dimnames(xx)[[2]]
     for (i in 1:nx) {
         for (j in 1:nx) {
             if (j < i) r[j, i] <- var(log(x[, i]/x[, j]))
             if (j == i) r[i, j] <- NA
             if (j > i) r[j, i] <- mean(log(x[, i]/x[, j]))
         }
     }
     r <- round(r, 3)
     for(i in 1:temp.x$m) r[i,i] <- NA
     cat(paste("Aitchison Variability Matrix, upper triangle contains variances ",
        "\nof log-ratios and lower triangle contains means of log-ratios,\nfor matrix ",
         deparse(substitute(xx)),", N = ", temp.x$n, "\n\n", sep = ""))
     print(r, na.print = " ")
     cat("\n")
     invisible()
}
