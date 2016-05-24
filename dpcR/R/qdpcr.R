#' Class \code{"qdpcr"}
#' 
#' An object representing digital PCR reaction depicted as Poisson process. Inherits 
#' from \code{\linkS4class{dpcr}}.
#' 
#' 
#' @name qdpcr-class
#' @aliases qdpcr-class qdpcr show.qdpcr show,qdpcr-method summary.qdpcr
#' summary,qdpcr-method
#' @docType class
#' @section Slots: \describe{ 
#' \item{list("mu")}{\code{"numeric"} of the expected number of events in a
#' defined interval.}
#' \item{list("C0")}{\code{"numeric"} the occurrence in a defined interval.}
#' \item{list("CT")}{\code{"numeric"} value of the
#' "average time" between the occurrence of a positive reaction and another
#' positive reaction.}
#' }
#' @author Stefan Roediger, Michal Burdukiewicz.
#' @seealso \code{\link{plot.qdpcr}},
#' @keywords classes
NULL

setClass("qdpcr", contains = "dpcr", representation(qpcr = "matrix", 
                                                    mu = "numeric", 
                                                    CO = "numeric",
                                                    CT = "numeric"))




#' Plot \code{qdpcr} objects
#' 
#' An analytical plot describing relationship between the cycle number and the
#' current value of Poisson mean. The plot can be used for quality control of
#' process.
#' 
#' The \code{rug} parameter allows user to add density of the number of events
#' to the plot.
#' 
#' @name plot.qdpcr
#' @aliases plot.qdpcr plot,qdpcr-method plot,qdpcr,ANY-method
#' @param x is a \code{\linkS4class{qdpcr}} object.
#' @param mincyc is the first cycle to start the plot from.
#' @param maxcyc the last cycle for the plot.
#' @param rug Adds a rug representation of the data to the plot.
#' @param digits how many significant digits are to be used in plot.
#' @author Stefan Roediger, Michal Burdukiewicz
#' @seealso \code{\linkS4class{qdpcr}}
#' @keywords hplot
#' @examples
#' 
#' library(qpcR)
#' test <- cbind(reps[1L:45, ], reps2[1L:45, 2L:ncol(reps2)], reps3[1L:45, 
#'         2L:ncol(reps3)])
#'         
#' plot(qpcr2pp(data = test, cyc = 1, fluo = NULL, model = l5, delta = 5), rug = TRUE)
#' 
NULL

setMethod("plot", signature(x = "qdpcr"), function(x, mincyc = 1, maxcyc = 45, rug = TRUE,
                                                   digits = getOption("digits") - 3) {
  # Plot the calculated qPCR data as Poisson processes
  res_qPCR <- slot(x, "qpcr")
  plot(res_qPCR[, 1], res_qPCR[,3], xlim = c(mincyc, maxcyc), 
       ylim = c(0, nrow(res_qPCR)), xlab = "Cycle", 
       ylab = expression(paste(lambda,
                               " (cycles)")), type = "S", lwd = 1.5)
  abline(h = nrow(res_qPCR) * 0.5, col = "grey")
  legend_texts <- c(paste0("n: ", slot(x, "n")),
                    paste0("k: ", sum(slot(x, "qpcr")[, "result"])),
                    as.expression(bquote(paste(mu, ": ", .(slot(x, "mu"))))),
                    paste0("CT: ", format(slot(x, "CT"), digits = digits)),
                    paste0("CO: ", format(slot(x, "CO"), digits = digits)))
  legend(mincyc, nrow(res_qPCR), legend_texts)
  # Add rug to the plot to illustrate the density of events
  if (rug) 
    rug(res_qPCR[, 1])
})



# setMethod("show", signature(object = "qdpcr"), function(object) {
#   print(slot(object, ".Data"))    
# })


# setMethod("summary", signature(object = "qdpcr"), function(object, print = TRUE) {
#   cat("\nmu: ", slot(object, "mu"), "\n")
#   cat("C0: ", format(slot(object, "CO")), "\n")
#   cat("Cycle time: ", format(slot(object, "CT")), "\n")
#   cat("Number of partitions: ", slot(object, "n"), "\n")
#   cat("Number of events: ", slot(object, "events"), "\n")
#   cat("\n")
# })
