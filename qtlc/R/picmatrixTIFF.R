#' Converting TIFF to matrix.
#'
#' Mostly internal function used by \code{createTLC} function. Additional parameters from \code{createTLC} goes there (RGB, comb).
#' @param ff TIFF file
#' @param RGB  Boolean, TRUE - keeps Red, Green and Blue intensities as three matrices. FALSE - using \code{comb} to combine intensities.
#' @param comb Vector, combines intensities according to luma. A vector containing three values for R, G and B conversion.
#'
#' @return Returns combined intensities matrix, or separated R, G, B matrices.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' #Internal function used by createTLC(...)
#' }
#'
#' @export
#' @importFrom tiff readTIFF
#'
picmatrixTIFF <- function (ff, RGB=TRUE, comb=c(0.30, 0.59, 0.11)) {
	
	rmat <- tiff::readTIFF(ff, convert=T);
	rmat <- 1 - rmat;
	
	#if pic is already GreyScale, then just returns the matrix
	if (dim(rmat)[3] == 1) {return(t(rmat));}
	
	#converts RGB pic to three matrices R,G, and B
	else if (RGB == TRUE) {
		
		red   <- t(rmat[,,1]);
		green <- t(rmat[,,2]);
		blue  <- t(rmat[,,3]);
		mat <- array(c(red, green, blue), dim=c(nrow(red), ncol(red), 3));
		cat("\nThe picture is RGB. Three matrices are created.\n");
		return(mat);
		}
		
		#converts to GreyScale with channels combining
		else if (RGB == FALSE) {
			red   <- rmat[,,1] * comb[1];
			green <- rmat[,,2] * comb[2];
			blue  <- rmat[,,3] * comb[3];
			
			mat <- red + green + blue;
			return(t(mat));
			}
	}

