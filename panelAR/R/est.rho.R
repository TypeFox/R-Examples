### Function to estimate rho
### Author: Konstantin Kashin
### August 1, 2013

# input: vector of residuals for a given panel
# k is rank
# output: correlation coefficient
# function is robust to panels with 1 time period
est.rho <- function(e,k,rhotype){
	mat <- embed(e,2)
	
	if(rhotype %in% c("breg","freg")){
		if(rhotype=="breg"){
		# backward regression
		# col 1 is resid, col 2 is lag
		num <- sum(apply(mat,MARGIN=1,prod),na.rm=TRUE)
		denom <- sum((mat[!is.na(mat[,1]),2])^2,na.rm=TRUE)
		rho <- num/denom
		} else {
		#forward regression
		# col 2 is resid, col 1 is forward
		num <- sum(apply(mat,MARGIN=1,prod),na.rm=TRUE)
		denom <- sum((mat[!is.na(mat[,2]),1])^2,na.rm=TRUE)
		rho <- num/denom
		}
	} else if(rhotype %in% c("dw","theil-nagar")){
		sse <- sum(e^2,na.rm=TRUE)
		sseN <- length(na.omit(e))
		dwal <- sum((mat[,1]-mat[,2])^2, na.rm=TRUE)/sse
		
		if(rhotype=="dw"){
		# Durbin-Watson calculation 
		rho <- 1-dwal/2
		} else{
		# Theil-Nagar
		rho <- (sseN^2 * (1-dwal/2) + k^2)/(sseN^2 - k^2)
		}	
	} else{
		sse <- sum(e^2,na.rm=TRUE)
		sseN <- length(na.omit(e))
		cov <- sum(apply(mat,MARGIN=1,prod),na.rm=TRUE)
		# scorr
		rho <- cov/sse
		if(rhotype=="theil"){
			# Theil rho
			# scale scorr by (T-k)/(T-1)
			# Ncov is T
			# Ncov is sseN-1 = T-1
			Ncov <- length(na.omit(apply(mat,MARGIN=1,prod)))
			rho <- (sseN - k)/Ncov*rho
		}
	}
	rho
}