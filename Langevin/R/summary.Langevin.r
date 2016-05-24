#' Summarize estimated drift and diffusion coefficients
#'
#' summary method for class "Langevin".
#'
#'
#' @param object an object of class "Langevin".
#' @param ... arguments to be passed to or from other methods. They are
#' ignored in this function.
#' @param digits integer, used for number formatting with \code{\link{signif}()}.
#'
#' @return The function \code{summary.Langevin()} returns a sumamry of the
#' estimated drift and diffusion coefficients
#'
#' @author Philip Rinn
#' @importFrom stats median
#' @export
summary.Langevin <- function(object, ..., digits=max(3, getOption("digits") - 3)) {
    cat(paste0(" Number of bins: ", dim(object$D1)[1],
               if(dim(object$D1)[2] > 1) paste0("x", dim(object$D1)[2]), "\n"),
        paste0("Population of the bins:\n",
              "\tMin.  : ", min(object$density, na.rm=T), "\n",
              "\tMedian: ", round(median(object$density, na.rm=T)), "\n",
              "\tMean  : ", round(mean(object$density, na.rm=T)), "\n",
              "\tMax.  : ", max(object$density, na.rm=T), "\n"),
        paste0("Number of NA's for D1: ", length(which(is.na(object$D1))),
               "\n"),
        paste0("Number of NA's for D2: ", length(which(is.na(object$D2))),
               "\n"),
        if(dim(object$D1)[2] == 1) paste0("Ratio between D4 and D2^2:\n",
              "\tMin.  : ", signif(min(object$D4/object$D2^2, na.rm=T), digits), "\n",
              "\tMedian: ", signif(median(object$D4/object$D2^2, na.rm=T), digits), "\n",
              "\tMean  : ", signif(mean(object$D4/object$D2^2, na.rm=T), digits), "\n",
              "\tMax.  : ", signif(max(object$D4/object$D2^2, na.rm=T), digits), "\n")
        )
}
