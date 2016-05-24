moment.gb2 <- function(k,shape1,scale,shape2,shape3){
	pk <- shape2+ k/shape1
	qk <- shape3- k/shape1
	if (qk <0) {print("moment does not exist: k >= aq", quote=FALSE);return(NA)}
	if (pk <0) {print("moment does not exist: k <= -ap", quote=FALSE);return(NA)}	
	logEk <- k*log(scale) + lgamma(pk) + lgamma(qk) - lgamma(shape2) - lgamma(shape3)
#	Ek <- (scale^k)*(gamma(pk)/gamma(shape2))*(gamma(qk)/gamma(shape3))
	return(exp(logEk))
}

incompl.gb2 <- function(x,k,shape1,scale,shape2,shape3) {
	pk <- shape2+ k/shape1
	qk <- shape3- k/shape1
	if (qk <0) {print("error: k >= aq", quote=FALSE);return(NA)}
	if (pk <0) {print("error: k <= -ap", quote=FALSE);return(NA)}
	return(pgb2(x,shape1,scale,pk,qk))
	}

el.gb2 <- function(shape1,scale,shape2,shape3) {log(scale)+(digamma(shape2)-digamma(shape3))/shape1}
vl.gb2 <- function(shape1,shape2,shape3) {(trigamma(shape2)+trigamma(shape3))/shape1^2}
sl.gb2 <- function(shape2,shape3) {(psigamma(shape2,deriv=2)-psigamma(shape3,deriv=2))/(vl.gb2(1,shape2,shape3))^(3/2)}
kl.gb2 <- function(shape2,shape3) {(psigamma(shape2,deriv=3)+psigamma(shape3,deriv=3))/(vl.gb2(1,shape2,shape3))^2}