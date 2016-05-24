gx.rqpca.print <-
function(save, ifload = TRUE, ifcntrb = FALSE, ifscore = TRUE)
{
     # Function to print the PCA matrices saved from gx.mva, gx.mva.closed,
     # gx.robmva, gx.robmva.closed and gx.rotate.  To save the matrices as
     # '.csv' files use function gx.rqpca.save.  The last table saved as an
     # object for easy access by users, rather than having to extract the
     # scores from a saved object.
     #
     cat("  PCA matrices for", deparse(substitute(save)), 
         "\n  Source data matrix:", save$input, "\n\n")
     if(ifload) {
         cat("  PC Loadings:\n")
         table.rows <- round(save$rload, 5)
         dimnames(table.rows)[[1]] <- save$matnames[[2]]
         dimnames(table.rows)[[2]] <- paste("PC-", 1:save$p, sep="")
         print(table.rows)
         cat("\n")
     }
     if(ifcntrb) {
         cat("  The relative % contributions of the scores of the loadings:\n")
         table.rows <- round(save$rcr, 2)
         dimnames(table.rows)[[1]] <- save$matnames[[2]]
         dimnames(table.rows)[[2]] <- paste("PC-", 1:save$p, sep="")
         print(table.rows)
         cat("\n  The cumulative relative % contributions of the PCs:\n")
         temp.table <- matrix(nrow=save$p, ncol=save$p)
         for (i in 1:save$p) {
             cumulation <- cumsum(save$rcr[i, ])
             temp.table[i, ] <- round(cumulation, 2)
         }
         dimnames(temp.table)[[1]] <- save$matnames[[2]]
         dimnames(temp.table)[[2]] <- paste("PC-", 1:save$p, sep="")
         print(temp.table)
         cat("\n")
     }
     if(ifscore) {
         cat("  Scores on the PCs:\n")
         table.rows <- round(save$rqscore, 4)
         dimnames(table.rows)[[2]] <- paste("PC-", 1:save$p, sep="")
         print(table.rows)
         cat("\n")
     }
     if(!is.null(save$nr)) {
         cat("  Following a Varimax rotation:\n\n")
         if(ifload) {
             cat("  Rotated PC Loadings:\n")
             table.rows <- round(unlist(unclass(save$vload)), 5)
             dimnames(table.rows)[[1]] <- save$matnames[[2]]
             dimnames(table.rows)[[2]] <- paste("RPC-", 1:save$nr, sep="")
             print(table.rows)
             cat("\n")
         }
         if(ifscore) {
             cat("  Scores on the Rotated PCs:\n\n")
             table.rows <- round(save$vscore, 4)
             dimnames(table.rows)[[2]] <- paste("RPC-", 1:save$nr, sep="")
             print(table.rows)
             cat("\n")
         }
     }
     invisible(table.rows)
}
