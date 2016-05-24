#' Retention factor (Rf)
#'
#' Calculates Rf values of the spots based on the marked start and stop of the solvent path.
#' @param object S3 object of the working TLC
#' @param sf Boolean, default FALSE - Start and Front should be marked. If TRUE, Start and Front were defined.
#'
#' @return Returns S3 object with new variables. \item{object$Rf_start}{Location of the solvent start on the TLC plate} \item{object$Rf_front}{Location of the solvent end on the TLC plate} \item{object$Rf}{Rf values of the spots}
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' #for more interactive variant; locate using mouse
#' object <- Rf(object, sf=TRUE)
#' }
#'
#' @export
#' @importFrom graphics locator abline rect text points strheight strwidth
#' 
Rf <- function(object, sf=F) { #sf - start, front F- false for auto

        mat <- object$mat;
        lok <- object$spots;
	st <- 0; fr <- ncol(mat); #automatic start/front
        
	if (sf == TRUE) {
		cat("\nLocate START, then FRONT.\n");
		l <- graphics::locator(type="p", pch=4, col="white", n=2);
		st <- l$y[1];
		fr <- l$y[2];
		graphics::abline(st, 0, col="red");
		graphics::rect(nrow(mat)/2-graphics::strwidth("START")/2, st+graphics::strheight("START")/2, nrow(mat)/2+graphics::strwidth("START")/2, st-graphics::strheight("START")/2, col="red", border=F);
		graphics::text(nrow(mat)/2, st, "START", col="white", cex=0.8);
		graphics::abline(fr, 0, col="red");
		graphics::rect(nrow(mat)/2-graphics::strwidth("START")/2, fr+graphics::strheight("START")/2, nrow(mat)/2+graphics::strwidth("START")/2, fr-graphics::strheight("START")/2, col="red", border=F);
		graphics::text(nrow(mat)/2, fr, "FRONT", col="white", cex=0.8);
		}
	
	RF <- c();	
	for (i in lok$y) {
		if (sf == TRUE) {poz <- i - st;} #izmena 28.02.2015.
			else {poz <- i};
		RF <- c(RF, poz/(fr-st));
            }
        object$Rf_start <- st;
        object$Rf_front <- fr;
	object$Rf <- RF;
	return(object);
}
