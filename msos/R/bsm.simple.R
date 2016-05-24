bsm.simple <-
function(x,y,z) {
	yz <- y%*%z%*%solve(t(z)%*%z)
	cx <- solve(t(x)%*%x)
	beta <- cx%*%t(x)%*%yz
	residz <- yz - x%*%beta
	df <- nrow(yz)-ncol(x)
	sigmaz <- t(residz)%*%residz/df
	se <- sqrt(outer(diag(cx),diag(sigmaz),"*"))
	list(Beta = beta, SE = se, T = beta/se, Covbeta = kronecker(cx,sigmaz),df = df,
		Sigmaz = sigmaz, Cx=cx)
}
