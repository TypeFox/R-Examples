#' Shows areadens2D plus background noise as segmented line
#'
#' Plots area-denses of the spot and background.
#' @param object S3 object of working TLC
#' @param spot Number of the spot (counted left to right).
#' @param plot Boolean, TRUE default and displays densitometric distribution.
#' @param returndf Boolean, TRUE by default, returns \code{data.frame} with \code{x} and \code{Area} values.
#' @param ptype Point type for the plot. Default "o" (Uses same values as \code{type} variable from \code{plot} function)
#' @param ... Additional parameters (for \code{plot} type function).
#'
#' @return Returns \code{data.frame} with \code{x} and \code{Area} values.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' #see areadens2D first
#' areadensnoisecompare2D(object, spot=3, ptype="l")
#' }
#'
#' @export
#' @importFrom graphics plot points
#' 
areadensnoisecompare2D <- function(object, spot=NULL, plot=TRUE, returndf=TRUE, ptype="o", ...) {

    if (is.null(spot)) {
        cat("\nSpecify the spot for density plot.\n");
        return(NULL);
    }
    if(is.null(object$noise_mat)) {
        cat("\nUsing noisepoly2D function is required before comparation!\n");
        return(FALSE);
    }

    mat <- object$spot_matrices[,,spot];
    nmat <- object$noise_mat[,,spot];
    P_original <- NULL;
    P_noiserem <- NULL;
    
    for(i in 1:nrow(mat)) {
        P_original <- c(P_original, sum(mat[i,]));
        P_noiserem <- c(P_noiserem, sum(nmat[i,]));
    }

    Pdf <- data.frame(1:nrow(mat), P_original, P_noiserem);
    names(Pdf) <- c("x", "AreaOriginal", "AreaNoiseRem");
    if (plot == TRUE) {
        #par(mfrow=c(1,2));
        graphics::plot(Pdf$x, Pdf$AreaOriginal, type=ptype, main=spot, xlab="x", ylab="Area", ...);
        graphics::points(Pdf$x, Pdf$AreaNoiseRem, type=ptype, lty=2, ...);
    }
    if (returndf == TRUE) {return(Pdf);}
}

