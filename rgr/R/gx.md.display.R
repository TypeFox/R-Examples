gx.md.display <-
function (xx, pcut = 0.1, ifprint = TRUE, file = NULL) 
{
     # Function to display Mahalanobis Distances and membership probabilities
     # together with selected variables from the dataframe or matrix used for
     # Mahalanobis distance estimations.  The data are sorted in order of
     # increasing probability of group membership and only those 'samples'
     # probabilities less that pcut are displayed.
     # Alternately, the entire table may be exported as a .csv file for later
     # use or display with a spreadsheet program.  If file is set to "" or
     # " " a default file name is generated.
     #
     # The dataframe from which the matrix passed for Mahalanobis Distance 
     # estimation was generated must be attached so that the data for the
     # variables to be appended to the Mahalanobis Distances and
     # probabilities of group membership are available.  In creating the
     # data matrix, xx, to be passed to the function with cbind the 'MD_s' 
     # and 'ppm_s' from the saved object MUST be in positions 1 and 2 for
     # the function to sort and display correctly.
     #
     dimnames(xx)[[2]][1:2] <- c("MD", "p_gm") 
     ppm <- xx[, 2]; nrows <- length(ppm[ppm < pcut])
     table.rows <- gx.sort(xx, 1, reverse = TRUE)
     table.rows[, 1:2] <- signif(table.rows[, 1:2], 3)
     #
     if(ifprint) {
         cat(paste("\n  Table of Mahalanobis Distances where probabilities of ", 
             "group membership (p_gm) are <", pcut, sep = ""), "\n\n")
         print(table.rows[1:nrows, ], print.gap = 2)
         cat("\n")
     }
     #
     if(!is.null(file)) {
         if(file == "" | file == " ") folder <- getwd()
         else folder <- file 
         filename <- paste(folder, "/MD_display.csv", sep="")
         write.csv(xx, file = filename, row.names = TRUE)
         cat("  Saved table will be in:\n  ", filename, "\n")
     }
     # 
     invisible(table.rows)
}
