
plotC <- function (object, columns = NULL) 
{
     if (is.null(object[["C"]])) stop("No random parameters to plot.")
  
     C <- object[["C"]]
     gNIV <- ncol(C) - 2
     numIDs <- nrow(C)
     graphics.off()
     par(ask = TRUE)
     
     if (is.null(columns)) columns <- 1:gNIV + 2
     
     for (i in columns) {
          plot(density(C[, i]), type = "l", main = paste(colnames(C)[i], "\n",
               sum(C[,i] >= 0), " (", round(sum(C[,i] >= 0) / length(C[,i]) * 100, 1),
               "%) >= 0", sep = ""), 
               xlab = "Utility", ylab = "Density")
               dev.flush()
     }
     par(ask = FALSE)
}