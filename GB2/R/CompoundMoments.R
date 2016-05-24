# Moment of order k
mkl.cgb2 <- function(k, shape1, scale, shape2, shape3, pl0, decomp="r"){
    if (decomp=="r") 
       {sh <- shape3
        a0 <- shape1}
    if (decomp=="l")
       {sh <- shape2
        a0 <- -shape1}
	Egb2 <- moment.gb2(k,shape1,scale,shape2,shape3)
	u2 <- qgamma(cumsum(pl0),sh)
	u1 <- c(0,u2[-length(pl0)])
	shk <- sh - k/a0  
    fac <- (pgamma(u2,shk) - pgamma(u1,shk))/(pgamma(u2,sh)-pgamma(u1,sh))
return(Egb2*fac)
}

# Moment of order k, -ap < k < aq
moment.cgb2 <- function(k, shape1, scale, shape2, shape3, pl0, pl, decomp="r"){
	pk <- shape2 + k/shape1
	qk <- shape3 - k/shape1
	if (qk <0) {print("moment does not exist: k >= aq", quote=FALSE);return(NA)}
	if (pk <0) {print("moment does not exist: k <= -ap", quote=FALSE);return(NA)}	
	Ek <- mkl.cgb2(k,shape1,scale,shape2,shape3,pl0,decomp)
	return(sum(pl*Ek))
}

# Incomplete moment of order k, -ap < k < aq
incompl.cgb2 <- function(x, k, shape1, scale, shape2, shape3, pl0, pl, decomp="r"){
	pk <- shape2+ k/shape1
	qk <- shape3- k/shape1
	if (qk <0) {print("moment does not exist: k >= aq", quote=FALSE);return(NA)}
	if (pk <0) {print("moment does not exist: k <= -ap", quote=FALSE);return(NA)}	
	    if (decomp=="r") 
         {sh <- shape3
          a0 <- shape1}
        if (decomp=="l")
         {sh <- shape2
          a0 <- -shape1}
        shk <- sh -k/a0
	u2 <- qgamma(cumsum(pl0),sh)
	ppl0 <- pgamma(u2,shk)
	ppl0 <- c(ppl0[1],diff(ppl0))
	Fk <- pl.cgb2(x,shape1,scale,pk,qk,ppl0,decomp) 
	Ek <- mkl.cgb2(k,shape1,scale,shape2,shape3,pl0,decomp)
	Mk <- Ek*Fk
	num <- sum(pl*Mk)
	denom <- sum(pl*Ek)
	return(num/denom)
}

