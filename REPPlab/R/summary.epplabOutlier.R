#' Summarize an epplabOutlier Object
#' 
#' Summarizes and prints an \code{epplabOutlier} object in an informative way.
#' 
#' The main information provided here is a table with names of the observations
#' which are considered outliers and in how many PP directions they are
#' considered outliers. This function is useful if the data has been given row
#' names.
#' 
#' @name summary.epplabOutlier
#' @aliases summary.epplabOutlier summary,epplabOutlier-method
#' @docType methods
#' @param object Object of class \code{epplabOutlier}.
#' @param ... Additional parameters
#' @author Klaus Nordhausen
#' @keywords methods print
#' @examples
#' 
#' # creating data with 3 outliers
#' n <-300 
#' p <- 10
#' X <- matrix(rnorm(n*p),ncol=p)
#' X[1,1] <- 9
#' X[2,4] <- 7 
#' X[3,6] <- 8
#' # giving the data rownames, obs.1, obs.2 and obs.3 are the outliers.
#' rownames(X) <- paste("obs",1:n,sep=".")
#' 
#' PP<-EPPlab(X,PPalg="PSO",PPindex="KurtosisMax",n.simu=20, maxiter=20)
#' OUT<-EPPlabOutlier(PP, k = 3, location = median, scale = mad)
#' summary(OUT)
#' 
#' @export
summary.epplabOutlier <- function(object, ...)
    {
    RowSUMS <- rowSums(object$outlier)
    RowMEANS <- RowSUMS/NCOL(object$outlier)*100
    
    if (sum(RowSUMS) !=0){
        IDrowSUMS <- RowSUMS>0
        OutlierNames <- row.names(object$outlier)[IDrowSUMS]
        OutlierFreq <- RowSUMS[IDrowSUMS]
        OutlierRelFreq <- RowMEANS[IDrowSUMS]
        OutlierN <- length(OutlierNames)
    
        OutlierTable <- as.table(rbind(OutlierNames, OutlierFreq, OutlierRelFreq))
        attr(OutlierTable,"dimnames")<-list(c("OutlierID: ", "Frequency: ", "Percentage:"),rep("",OutlierN))

        cat("REPPlab Outlier Summary\n")
        cat("-----------------------\n")
        cat("Index name       :", object$PPindex, "\n")
        cat("Algorithm used   :", object$PPalg, "\n")
        cat("Location used    :", object$location, "\n")
        cat("Scale used       :", object$scale, "\n")
        cat("k value used     :", object$k, "\n")
        cat("-----------------------\n")
        cat("\n")
        cat("Number of outliers detected:\n", length(OutlierNames), "\n")
        cat("\n")
        cat("Observations considered outliers:")
        print(OutlierTable)
        cat("\n")} else {
        cat("REPPlab Outlier Summary\n")
        cat("-----------------------\n")
        cat("Index name       :", object$PPindex, "\n")
        cat("Algorithm used   :", object$PPalg, "\n")
        cat("Location used    :", object$location, "\n")
        cat("Scale used       :", object$scale, "\n")
        cat("k value used     :", object$k, "\n")
        cat("-----------------------\n")
        cat("\n")
        cat("No outliers detected")
        cat("\n")
        }
    invisible(object)
    }
