writeModel <- function(object, writeDraws = FALSE, path = getwd()) {
     
     # Error check
     if (!"RSGHB" %in% class(object)) stop("'object' is not of class RSGHB")
     
     # Disk location
     orig.path <- getwd()
     setwd(path)
     
     gSIGDIG <- object[["gSIGDIG"]]
     modelname <- object[["modelname"]]
     
     # Make sure the output files don't already exist
     # If they do, append ~+1 to the model name.
     orig <- modelname
     i <- 1
     # If we add more file types to this function, need to add them here
     while (any(file.exists(paste0(modelname, c(".log", "_A.csv", "_B.csv", "_Bsd.csv", "_C.csv", "_Csd.csv", "_D.csv", "_F.csv")))))
     {
          modelname <- paste0(orig, "~", i)
          i <- i + 1
     }
     
     if (modelname != orig) warning(paste0("Model files associated with '", orig, "' already exist. Writing results as '", modelname, "'"))
     
     ### Write out files for this model object if available
     orig.options <- options() # Set print width wide enough for the log file
     options(width = 1000)
     
     # log file
     sink(paste0(object[["modelname"]], ".log"))
     cat("Model Name:", object[["modelname"]], "\n", sep = "\t")
     cat("Number of individuals:", object[["gNP"]], "\n", sep = "\t")
     cat("Number of observations:", object[["gNOBS"]], "\n", sep = "\t")
     cat("Number of preliminary iterations:", object[["gNCREP"]], "\n", sep = "\t")
     cat("Number of draws used per individual:", object[["gNEREP"]], "\n", sep = "\t")
     cat("Random Seed:", object[["gSeed"]], "\n", sep = "\t")
     cat("Total iterations:", object[["gNCREP"]] + object[["gNEREP"]], "\n", sep = "\t")
     
     if (!is.null(object[["df"]])) cat("Degrees of Freedom:", object[["df"]], "\n", sep = "\t")
     
     cat("Number of parameters:", length(object$params.vary) + length(object$params.fixed), "\n\n", sep = "\t")
     
     if (length(object[["params.fixed"]]) > 0)
     {
          cat("Fixed parameters estimated:\n")
          cat(paste0(object[["params.fixed"]], "\n", collapse = ""))
     }
     cat("\n")
     if (length(object[["params.vary"]]) > 0)
     {
          cat("Random parameters estimated (Distribution):", "\n")
          cat(paste0(paste(object[["params.vary"]], "(", object[["distributions"]], ")"), "\n"), collapse = "", sep = "")
     }
     cat("\n")
     if (!is.null(object[["constraints"]]))
     {
          cond <- c("<", ">")
          cat("Constraints applied to random parameters (param1 - inequality - param2):\n")
          for(i in 1:length(object[["constraints"]]))
          {
               if(object[["constraints"]][[i]][3] == 0) cat(object[["params.vary"]][i], cond[object[["constraints"]][[i]][2]], 0, "\n")
               if(object[["constraints"]][[i]][3] != 0) cat(object[["params.vary"]][i], cond[object[["constraints"]][[i]][2]], object[["params.vary"]][object[["constraints"]][[i]][3]], "\n")
          }
          
     }
     cat("\n")
     if (!is.null(object[["pv"]])) {
          cat("Prior Variance-Covariance Matrix:\n")
          print(object[["pv"]])
     }
     
     cat("\n-----------------------------------------------------------\n\n")
     
     print(object[["iter.detail"]], row.names = FALSE)
     
     sink()
     
     options(orig.options)
     
     # prior variance matrix
     if (!is.null(object[["pv"]]))  write.table(object[["pv"]], paste0(modelname, "_pvMatrix.csv"), sep = ",", col.names = NA)
     
     # CSVs
     if (!is.null(object[["A"]]))   write.table(signif(object[["A"]],   gSIGDIG), paste0(modelname, "_A.csv"),   sep = ",", row.names = FALSE)
     if (!is.null(object[["B"]]))   write.table(signif(object[["B"]],   gSIGDIG), paste0(modelname, "_B.csv"),   sep = ",", row.names = FALSE)
     if (!is.null(object[["Bsd"]])) write.table(signif(object[["Bsd"]], gSIGDIG), paste0(modelname, "_Bsd.csv"), sep = ",", row.names = FALSE)
     if (!is.null(object[["C"]]))   write.table(signif(object[["C"]],   gSIGDIG), paste0(modelname, "_C.csv"),   sep = ",", row.names = FALSE)
     if (!is.null(object[["Csd"]])) write.table(signif(object[["Csd"]], gSIGDIG), paste0(modelname, "_Csd.csv"), sep = ",", row.names = FALSE)
     if (!is.null(object[["F"]]))   write.table(signif(object[["F"]], gSIGDIG), paste0(modelname, "_F.csv"), sep = ",", row.names = FALSE)
     
     # D file
     if (!is.null(object[["D"]])) {
          
          # 2D labeling
          labelmatrix <- matrix(1:(length(object[["params.vary"]])^2), length(object[["params.vary"]]), length(object[["params.vary"]]))
          rownames(labelmatrix) <- colnames(labelmatrix) <- object[["params.vary"]]
          dlabels <- paste0(rownames(labelmatrix)[row(labelmatrix)[lower.tri(labelmatrix, diag = TRUE)]],
                            " x ",
                            colnames(labelmatrix)[col(labelmatrix)[lower.tri(labelmatrix, diag = TRUE)]])
          md.write <- matrix(0, nrow = dim(object[["D"]])[3], ncol = length(dlabels) + 1)
          colnames(md.write) <- c("Iteration", dlabels)
          
          # Convert to 2D
          md.write[, "Iteration"] <- as.numeric(dimnames(object[["D"]])[[3]])
          for (i in 1:dim(object[["D"]])[3]) md.write[i, -1] <- object[["D"]][,,i][lower.tri(object[["D"]][,,i], diag = TRUE)]
          
          write.table(signif(md.write, gSIGDIG), paste0(modelname, "_D.csv"), sep = ",", row.names = FALSE)
     }
     
     # individual draws
     if (!is.null(object[["Draws"]]) & writeDraws) {
          cat("Creating individual draw files, this may take a few minutes.\n")     
          for(i in 1:length(object[["Draws"]]))
          {
               
               fn <- paste0("Draws_", names(object[["Draws"]])[[i]], ".csv") 
               write.table(signif(object[["Draws"]][[i]], gSIGDIG), fn, sep = ",", row.names = FALSE, col.names = TRUE)
               
          }
     }
     
     setwd(orig.path)

}
