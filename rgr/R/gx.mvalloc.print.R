gx.mvalloc.print <-
function(save, ifprint = TRUE, unalloc = TRUE, file = NULL)
{
     # Function to print posterior allocations and predicted probabilities
     # of group membership from an object saved from gx.mvalloc or
     # gx.mvalloc.closed either to the screen or a file in the Working
     # Directory, or elsewhere, e.g.,'d:\\stuff\\ecg\\'.
     #
     kk <- save$kk
     n <- save$n
     cat("  Reference Groups:\n")
     for(k in 1:kk) cat("    ", k, "  ", save$groups[k], "\n")
     cat("  pcrit was set to:", save$pcrit, "\n\n")
     pgm <- save$pgm
     ix <- seq(1:n)
     if(ifprint) {
         cat("  Row    Group  Probabilities of Group Membership:\n\n")
         for(i in 1:n) {
             if(unalloc & save$xalloc[i] == 0) cat("  ", ix[i], "    ", 
                 save$xalloc[i], "  ", pgm[i, ], "\n")
             if(!unalloc) cat("  ", ix[i], "    ", save$xalloc[i], "  ", 
                 pgm[i, ], "\n")
         }
         cat("\n")
     }
     if(!is.null(file)) {
         if(file == "" | file == " ") folder <- getwd()
         else folder <- file
         save.name <- deparse(substitute(save))
         filename <- paste(folder, "/", save.name, "_mvalloc.csv", sep = "") 
         cat("  Saved table will be in:\n  ", filename, "\n\n")
         sink(filename)
         on.exit(sink())
         cat("Row,Group", paste(",pgm[", 1:kk, "]", sep = ""), sep = "")
         for(i in 1:n) cat("\n", i, ",", save$xalloc[i], paste(",", pgm[i, ], sep = ""), sep = "")
     }
     invisible()
}
