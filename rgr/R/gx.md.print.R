gx.md.print <-
function(save, pcut = 0.1, ifprint = TRUE, file = NULL) 
{
     # NOTE: When printed the Mahanalobis Distances are listed in descending
     # order, i.e. most extreme first.
     #
     save.name <- deparse(substitute(save))
     cat("  Mahalanobis Distances for", save.name, "\n  Source data matrix:",
         save$input, "\n\n")
     ppm <- save$ppm; nrows <- length(ppm[ppm < pcut])
     md <- save$md; row.i <- rev(order(md))
     table.rows <- cbind(save$matnames[[1]], signif(md, 3), signif(ppm, 3))
     dimnames(table.rows)[[2]] <- c("Row ID", "MD", "p_gm")
     #
     if(ifprint) {
         cat(paste("  Table of Mahalanobis Distances where probabilities of ", 
             "group membership (p_gm) are <", pcut, sep = ""), "\n\n",
             "  ID\t\t  MD     p_gm\n\n")
         for (i in 1:nrows) {
             cat("  ", table.rows[row.i[i], 1], "\t\t", 
                 table.rows[row.i[i], 2], "\t", table.rows[row.i[i], 3], "\n")
         }
         cat("\n")
     }
     if(!is.null(file)) {
         if(file == "" | file == " ") folder <- getwd()
         else folder <- file 
         filename <- paste(folder, "/", save.name, "_MDs.csv", sep="")
         cat("  Saved table will be in:\n  ", filename, "\n\n")
         write.csv(table.rows, file = filename, row.names = FALSE)
     }
     invisible(table.rows)
}
