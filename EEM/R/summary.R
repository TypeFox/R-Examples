#' SummarizeEEM EEM list
#' 
#' Summarize by listing the sample number, names and their dimensions
#' 
#' @param object a list containing EEM data as created by \code{readEEM} function.
#' @param ... arguments for \code{summary} function
#' 
#' @return Text on console
#' 
#' @examples
#' data(applejuice)
#' summary(applejuice)
#' 
#' @export

summary.EEM <-
    function(object, ...){
        
        # number of samples
        writeLines(paste("Number of samples:", length(object)))
        
        # sample name
        writeLines("Sample names: ")
        print(names(object))
        
        # dimension
        dimension <- sapply(object, dim)
        dimensionVector <- apply(dimension, 2, paste, collapse = "x")
        dimensionUnique <- unique(dimensionVector)
        dimN <- length(dimensionUnique)
        if (dimN == 1) {
            writeLines(paste("Dimension [EmxEx]:", dimensionUnique))
        } else {
            writeLines("Dimension [EmxEx]:")
            for (i in 1:dimN) {
                writeLines(paste("[", i, "] ", dimensionUnique[i], sep = ""))
            }
        }
        
        # EX range
        EX_range <- sapply(object, function(x) range(as.numeric(colnames(x))))
        unique_EX_range <- unique(t(EX_range))
        dim_EX_range <- nrow(unique_EX_range)
        if (dimN == 1) {
            writeLines(paste0("EX range: ", unique_EX_range[1], "~", unique_EX_range[2], " [nm]"))
        } else {
            writeLines("EX range: ")
            for (i in 1:dim_EX_range) {
                writeLines(paste0("[", i, "] ", unique_EX_range[i, 1], "~", unique_EX_range[i, 2]))
            }
        }
        
        # EM range
        EM_range <- sapply(object, function(x) range(as.numeric(rownames(x))))
        unique_EM_range <- unique(t(EM_range))
        dim_EM_range <- nrow(unique_EM_range)
        if (dimN == 1) {
            writeLines(paste0("EM range: ", unique_EM_range[1], "~", unique_EM_range[2], " [nm]"))
        } else {
            writeLines("EM range: ")
            for (i in 1:dim_EM_range) {
                writeLines(paste0("[", i, "] ", unique_EM_range[i, 1], "~", unique_EM_range[i, 2]))
            }
        }
        
    }
