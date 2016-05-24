`coef.mlcm` <- function(object, ...){
	if (object$method == "glm")
		cc <- object$obj$coef else {
		cc <- object$pscale[-c(1, length(object$stimulus) + 
			1)]
		dn <- as.vector(sapply(names(object$data)[c(2, 4)], 
				substr, 1, 1))
		if (object$model == "add")
		names(cc) <- as.vector(t(outer(dn, 2:((length(cc) + 2)/2), 
				paste, sep = ""))) else
		names(cc) <- paste(dn[object$whichdim], 2:(length(cc) + 1), 
			sep = "")
				}
		cc
}