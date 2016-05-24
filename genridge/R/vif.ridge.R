# variance inflation factors for ridge regression objects

vif.ridge <- function(mod, ...) {

	Vif <- function(v) {
		R <- cov2cor(v)
		detR <- det(R)
		p <- nrow(v)
		res <- rep(0,p)
		for (i in 1:p) {
			res[i] <- R[i,i] * det(as.matrix(R[-i,-i])) / detR
		}
		res
	}

#	if(!require("car")) stop("Requires the car package for the vif generic")
	V <- vcov(mod)
	res <- t(sapply(V, Vif))
	colnames(res) <- colnames(coef(mod))
	rownames(res) <- rownames(coef(mod))
	res
}
