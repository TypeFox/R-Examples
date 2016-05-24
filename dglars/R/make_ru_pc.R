make_ru_pc <- function(object){
	ru <- object$ru
	np <- object$np
	p <- object$p
	nav <- object$nav
	A <- object$A[1:nav]
	ru <- matrix(ru[1:(p*np)],p,np)
	rownames(ru) <- colnames(object$X)
	ru[sort(A),]
}
