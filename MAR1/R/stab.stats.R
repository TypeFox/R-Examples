stab.stats<-function(B,sigma,V_inf){

## RESILIENCE ##

# variance of the stationary distribution (time spent away from mean)
	# (higher means less stable; species interactions greatly amplify environmental variance)

	# eigenvalues of the B matrix

		eig.b<-eigen(B)$values

	# determinant of the B matrix

		det.b<-abs(det(B))^(2/nrow(B))

# return rate of the transition distribution to the stationary distribution (time to return to mean)
	# (higher means less stable; slower return rate)

	# asymmptotic rate of return of the mean
	# max eigenvalue of B matrix
		max.eig<-suppressWarnings(max(as.numeric(eig.b)))

	# asymptotic rate of return of the variance
	# max eigenvalue of B matrix kronecker products

		max.eig.kr<-suppressWarnings(max(as.numeric(eigen(kronecker(B,B))$values)))

# REACTIVITY

# how much values move towards stationary distribution between time steps ("pull" towards mean)
# (higher means less stable; less tendency towards mean between time steps)

	# -tr(sigma)/tr(Vinf)

		covar_sigma_Vinf<--sum(diag(sigma))/sum(diag(V_inf))

	# max eigenvalue of B'B matrix ("worst-case" reactivity)

		max.eig.tbxb<-(max(eigen(t(B)%*%B)$values)-1)


list(
	resilience=list(
		eigB		=	eig.b,
		detB		=	det.b,
		maxeigB	=	max.eig,
		maxeigkrB	=	max.eig.kr),
	reactivity=list(
		sigma.over.Vinf = covar_sigma_Vinf,
		maxeigBxB	=	max.eig.tbxb)
)
		
# save all objects in run.mar function frame
# fun.obj<-ls()
# for(i in 1:length(fun.obj)){
# assign(fun.obj[i],get(fun.obj[i]),envir=parent.frame())}

}

