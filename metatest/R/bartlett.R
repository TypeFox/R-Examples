# bartlett <-
# function(w, x, partar, itar, n, p) {
# 	
# 	wsq=w^2
# 	imat=matrix(diag(c(rep(1,n))),n,n)                              
# 	xtar=matrix(x[1:n,itar],n,1)
# 	if(p>1)  r=x[,-itar]
# 	if(p==1) r=diag(c(rep(1,n))) # for intercept only model
# 	aa=solve(t(r)%*%w%*%r) 
# 	z=(imat-(r %*% aa %*% t(r) %*% w)) %*% xtar
# 	
# 	l1=-t(z) %*% w %*% z
# 	l2=-(0.5 * sum(wsq)+(partar^2) %*% t(z) %*% wsq %*% r %*% aa %*% t(r) %*% wsq %*% z)^(-1)
# 	l3=-2.0 * t(z) %*% w %*% (w - w %*% r %*% aa %*% t(r) %*% w) %*% w %*% z
# 	l4=t(z) %*% wsq %*% z
# 	
# 	scalfact = (1+0.5*(l1^(-1))*l2*(-l3 + 0.5*(l1^(-1))*(l4^2) + l4*(p-1)))^(-1)
# 	
# 	return(scalfact)
# 	
# }



bartlett <- function (w, x, partar, itar, n, p)
{
	wsq = w^2
	imat = matrix(diag(c(rep(1, n))), n, n)
	xtar = matrix(x[1:n, itar], n, 1)
	if (p > 1){
		r = x[, -itar]
		aa = solve(t(r) %*% w %*% r)
		z = (imat - (r %*% aa %*% t(r) %*% w)) %*% xtar
		l1 = -t(z) %*% w %*% z
		l2 = -(0.5 * sum(wsq) + (partar^2) %*% t(z) %*% wsq %*% r %*%
			aa %*% t(r) %*% wsq %*% z)^(-1)
		l3 = -2 * t(z) %*% w %*% (w - w %*% r %*% aa %*% t(r) %*%
			w) %*% w %*% z
		l4 = t(z) %*% wsq %*% z
		scalfact = (1 + 0.5 * (l1^(-1)) * l2 * (-l3 + 0.5 * (l1^(-1)) *
				(l4^2) + l4 * (p - 1)))^(-1)
	}
	if(p==1) {
		z=xtar
		l1=-t(z) %*% w %*% z
		l2 = -(0.5 * sum(wsq))^(-1)
		l3 = -2 * t(z) %*% w %*% (w ) %*% w %*% z
		l4 = t(z) %*% wsq %*% z
		scalfact = (1 + 0.5 * (l1^(-1)) * l2 * (-l3 + 0.5 * (l1^(-1)) *
				(l4^2)))^(-1)
	}
	return(scalfact)
}
