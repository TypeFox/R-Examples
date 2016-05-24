# 
# Ingmar Visser, oct 9, 2014
# 
# Compute a corrected hessian using
# 
hess2cov <- function(hess,lincon) {
	np <- dim(hess)[1]
	if(!(dim(hess)[1]==dim(hess)[2])) stop("'hess' should be a square matrix")
	se=rep(0,np)
	hs=hess	
	if(nrow(lincon)>0) {
		A=lincon
		d=hs+t(A)%*%A 
		di=try(solve(d),silent=TRUE)
		if(class(di)=="try-error") {
			warning("Hessian singular, ses could not be computed.") 
			val=list(se=0,hs=0)
		} else {
			ada=A%*%di%*%t(A)
			adai=try(solve(ada),silent=TRUE)
			if(class(adai)=="try-error") {
				warning("Near-singular hessian, ses may be bad.\n")
				diag(ada)=diag(ada)*1.000001
				adai=try(solve(ada))
				if(class(adai)=="try-error") {
					warning("Corrected hessian also singular, ses computed without contraints.\n")
					se=sqrt(diag(di))
				} else {
					ch=di-di%*%t(A)%*%adai%*%A%*%di
					se=sqrt(diag(ch))
				}
			} else {
				ch=di-di%*%t(A)%*%adai%*%A%*%di
				se=sqrt(diag(ch))
			} 
		} 
	} else {
		se=sqrt(diag(solve(hs)))
	}
	val=list(se=se,hs=hs)
	val
}