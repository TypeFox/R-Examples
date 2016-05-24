##########################################################
# general functions of H_I and Q_I  B&L Eqn 2.3, 2.4
# as in Table 1
# These are simply experimental, probably not to be exported
##########################################################


Jtr <- function (H, Q, a, b, f) {
		I <- diag(nrow(H))
		res <- H %*% Q %*% mpower(I-H-Q, a) %*% mpower(I-H, b)
		f * tr(res)
	}
Jdet <- function (H, Q, a, b, f) {
		I <- diag(nrow(H))
		res <- H %*% Q %*% mpower(I-H-Q, a) %*% mpower(I-H, b)
		f * det(res)
	}

	# Cook D, in terms of Jtr()
COOKD <- function(H, Q, n, p, r, m) {
 		f <- (n-p)/p
 		Jtr(H, Q, 0, -2, f)
 }

	# DFFITS^2, in terms of Jtr()
DFFITS <- function(H, Q, n, p, r, m) {
 		f <- (n-p-m)/p
 		Jtr(H, Q, -1, 0, f)
 }
 
COVRATIO <- function(H, Q, n, p, r, m) {
 		f <- ((n-p)/(n-p-m))^r*p
 		Jdet(H, Q, p, -(r+p), f)
 }
 
