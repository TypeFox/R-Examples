coef.haplin <- function(object, cov.type = "Fisher.EM", ...){
##
## EXTRACTS COEFFICIENTS AND THE VARIANCE-COVARIANCE MATRIX FROM A haplin OBJECT
##
#
## CHECK
if(!is.element(cov.type, c("Fisher.EM", "Fisher", "resamp", "robust"))) stop("Unknown covariance matrix requested!")
#
## EXTRACT GLM ESTIMATION RESULT
.res <- object$result
#
##
if(cov.type == "robust"){
	#
	warning("Computation of the robust covariance matrix has not been extensively tested!")
	## COMPUTE COEFFICIENTS, WITH "STANDARD" COVARIANCE MATRIX
	.ut <- coef.tri.glm(.res, cov.type = "Fisher.EM")
	###.score <- object$score$score # DENNE ER IKKE RIKTIG! REGNET UNDER NULL!
	.score <- object$temp1$score
	## COMPUTE ROBUST VARIANCE
	#.cov <- .ut$cov %*% t(.score) %*% .score %*% .ut$cov
	.cov <- crossprod(tcrossprod(.score, .ut$cov))
	## REPLACE SIMPLE WITH ROBUST
	attr(.cov, "cov.type") <- "robust"
	.ut$cov <- .cov
}else{
	.ut <- coef.tri.glm(.res, cov.type = cov.type)
}
#
##
return(.ut)
}
