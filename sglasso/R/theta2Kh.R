theta2Kh <- function(theta, mask, Dimnames){
	p <- dim(mask)[1]
	th_nm <- names(theta)
	Kh <- mask[upper.tri(mask, diag = TRUE)]
	id <- match(Kh, th_nm)
	Kh <- theta[id]
	if(any(is.na(Kh)))	Kh[is.na(Kh)] <- 0
	Kh <- new("dspMatrix", x = Kh, Dim = c(p, p), Dimnames = Dimnames)
	Kh
}