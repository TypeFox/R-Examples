gx.rqpca.save <-
function(save, ifload = TRUE, ifcntrb = FALSE, ifscore = TRUE, file = NULL)
{
     # Function to export the PCA matrices saved from gx.mva, gx.mva.closed,
     # gx.robmva, gx.robmva.closed and gx.rotate as '.csv' files to the 
     # Working Directory, the default, or to a folder provided by the user
     # in 'file'.  
     #
     if(is.null(file)) folder <- getwd()
     else folder <- file
     save.name <- deparse(substitute(save))
     cat("  PCA matrices for", save.name, "\n  Source data matrix:",
         save$input, "\n\n")
     #
     if(ifload) {
         table.rows <- round(save$rload, 5)
         dimnames(table.rows)[[1]] <- save$matnames[[2]]
         dimnames(table.rows)[[2]] <- paste("PC-", 1:save$p, sep="")
         filename <- paste(folder, "/", save.name, "_PC-loadings.csv",
             sep = "")
         write.csv(table.rows, file = filename, row.names = TRUE)
         cat("  PC loadings have been saved in:\n  ", filename, "\n\n")
     }
     if(ifcntrb) {
         table.rows <- round(save$rcr, 2)
         dimnames(table.rows)[[1]] <- save$matnames[[2]]
         dimnames(table.rows)[[2]] <- paste("PC-", 1:save$p, sep="")
         filename <- paste(folder, "/", save.name, "_PC-contributions.csv",
             sep = "")
         write.csv(table.rows, file = filename, row.names = TRUE)
         cat("  PC contributions have been saved in:\n  ", filename, "\n\n")
     }
     if(ifscore) {
         table.rows <- round(save$rqscore, 4)
         dimnames(table.rows)[[2]] <- paste("PC-", 1:save$p, sep="")
         filename <- paste(folder, "/", save.name, "_PC-scores.csv",
             sep = "")
         write.csv(table.rows, file = filename, row.names = TRUE)
         cat("  PC scores have been saved in:\n  ", filename, "\n\n")
     }
     if(!is.null(save$nr)) {
         if(ifload) {
             table.rows <- round(unlist(unclass(save$vload)), 5)
             dimnames(table.rows)[[1]] <- save$matnames[[2]]
             dimnames(table.rows)[[2]] <- paste("RPC-", 1:save$nr, sep="")
             filename <- paste(folder, "/", save.name,
                 "_Rotated-PC-loadings.csv", sep = "")
             write.csv(table.rows, file = filename, row.names = TRUE)
             cat("  Rotated PC loadings have been saved in:\n  ", filename,
                 "\n\n")
         }
         if(ifscore) {
             table.rows <- round(save$vscore, 4)
             dimnames(table.rows)[[2]] <- paste("RPC-", 1:save$nr, sep="")
             filename <- paste(folder, "/", save.name,
                 "_Rotated-PC-scores.csv", sep = "")
             write.csv(table.rows, file = filename, row.names = TRUE)
             cat("  Rotated PC scores have been saved in:\n  ", filename,
                 "\n\n")
         }
     }
     invisible(table.rows)
}
