#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
checkEmbParms <- function(series, m, d, t=0, s=1, ref=NULL) {
	n <- length(series)-(m-1)*d
	if(n <= 0) stop("Not enough points to handle these parameters")
	if(!is.null(ref)) if (ref>n) stop("Not enough points to handle these parameters")
	if(t<0) stop("Theiler window t must be non-negative")
	if(s<=0) stop("Number of steps must be positive")

}
