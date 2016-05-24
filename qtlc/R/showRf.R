#' Shows Rf on the plot
#'
#' Shows prior analysed Rf on the new plot of the 2D matrix.
#' @param object S3 object of the working TLC
#' @param col Color of the lines.
#' @param adjust Adjustment for the space of the text. Default value is usualy just OK.
#' @param cex A zoom factor for the text.
#'
#' @return None.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' showRf(object)
#' }
#'
#' @export
#' @importFrom graphics locator abline rect text points strheight strwidth
#' 
showRf <- function(object, col = "green", adjust = NULL, cex = 0.6) {

    mat <- object$mat;
    st <- object$Rf_start;
    fr <- object$Rf_front;
    val <- object$Rf;
    if (is.null(adjust)) { adjust <- round(0.065 * ncol(mat), 0);}

    graphics::abline(st, 0, col="red");
    graphics::rect(nrow(mat)/2-graphics::strwidth("START")/2, st+graphics::strheight("START")/2, nrow(mat)/2+graphics::strwidth("START")/2, st-graphics::strheight("START")/2, col="red", border=F);
    graphics::text(nrow(mat)/2, st, "START", col="white", cex=0.8);
    graphics::abline(fr, 0, col="red");
    graphics::rect(nrow(mat)/2-graphics::strwidth("START")/2, fr+graphics::strheight("START")/2, nrow(mat)/2+graphics::strwidth("START")/2, fr-graphics::strheight("START")/2, col="red", border=F);
    graphics::text(nrow(mat)/2, fr, "FRONT", col="white", cex=0.8);

    for(ts in 1:length(object$spots$y)) {
        graphics::points(c(object$spots$x[ts], object$spots$x[ts]), c(st, object$spots$y[ts]), type="l", col=col, lty=2);
    }
   
    graphics::text(object$spots$x, object$spots$y + adjust, round(object$Rf, 2), cex=cex, col="green");
}
