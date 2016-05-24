#' Polynomial estimation of the image noise.
#'
#' Estimation, and noise removal using polynomial model.
#' @param object S3 object of the working TLC
#' @param gd Defines position of the center of the rectangular samples of the image background (above or bellow located spots).
#' @param power Order of the polynome.
#' @param col Color of the borders of the rectangles for bkg samples.
#'
#' @return Returns S3 object with new variables. \item{object$noise_mat}{The 3D matrix (width, height, number of spots)} \item{object$noise_fit}{Linear model for the polynomial fit} \item{object$noisefit_spot_sums}{Sums of the noise samples areas}
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' #object is a tlc with 14 spots, and selection 50x80(wxh)
#' object <- noisepoly2D(object, gd=80, power=13)
#' }
#'
#' @export
#' @importFrom graphics rect plot
#' @importFrom grDevices x11
#' @importFrom stats lm predict poly
#' 
noisepoly2D <- function(object, gd=20, power=5, col="green") {
    
    lok <- object$spots;
    w <- object$mat_cell$w;
    h <- object$mat_cell$h;

    x1 <- round(lok$x-w/2);
    y1 <- round(gd + lok$y+h/2);
    x2 <- round(lok$x+w/2);
    y2 <- round(gd + lok$y-h/2);
    
    graphics::rect(x1, y1, x2, y2, lty=2, border=col);

    elements <- c();
    
    for(i in 1:length(lok$x)) {
        elements <- c(elements, object$mat[x1[i]:x2[i], y2[i]:y1[i]]);
    }
    
    noise_mat <- array(elements, dim=c(w, h, length(lok$x)));

    sums <- c();
    for(i in 1:dim(noise_mat)[3]) {
        sums <- c(sums, sum(noise_mat[,,i]));
    }

    noise_fit <- stats::lm(sums ~ stats::poly(lok$x, power, raw=TRUE));
    grDevices::X11();
    graphics::plot(lok$x, sums, main="Background noise polynomal fit", xlab="x", ylab="Signal");
    graphics::points(lok$x, stats::predict(noise_fit), type="o", pch=8, lwd=2, lty=2, col="red");

    nnoise <- c();
    for(jj in 1:length(stats::predict(noise_fit))) { nnoise <- c(nnoise, stats::predict(noise_fit)[[jj]]); }
    nnoise_spot_sums <- object$spot_sums - nnoise;

    object$noise_mat <- noise_mat;
    object$noise_fit <- noise_fit;
    object$noisefit_spot_sums <- nnoise_spot_sums;
    return(object);
}

