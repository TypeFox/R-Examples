#' Print estimated drift and diffusion coefficients
#'
#' print method for class "Langevin".
#'
#'
#' @param x an object of class "Langevin".
#' @param digits integer, used for number formatting with \code{\link{signif}()}.
#' @param ... further arguments to be passed to or from other methods. They are
#' ignored in this function.
#'
#' @return The function \code{print.Langevin()} returns an overview of the
#' estimated drift and diffusion coefficients.
#'
#' @author Philip Rinn
#' @export
print.Langevin <- function(x, digits=max(3, getOption("digits")-3), ...) {
    cat(paste0("Drift and diffusion estimates for the ",
               ifelse(dim(x$D1)[2] > 1, "two", "one"),
               " dimensional Langevin Approach\n"),
        paste0("Number of bins: ", dim(x$D1)[1],
               if(dim(x$D1)[2] > 1) paste0("x", dim(x$D1)[2]), "\n"),
        paste0("Range of the bins: ", signif(range(x$U,na.rm=T)[1], digits),
               " ... ", signif(range(x$U,na.rm=T)[2], digits),"\n"),
        paste0("Range of D1: ", signif(range(x$D1,na.rm=T)[1], digits),
               " ... ", signif(range(x$D1,na.rm=T)[2], digits),"\n"),
        paste0("Range of D2: ", signif(range(x$D2,na.rm=T)[1], digits),
               " ... ", signif(range(x$D2,na.rm=T)[2], digits),"\n")
        )
}
