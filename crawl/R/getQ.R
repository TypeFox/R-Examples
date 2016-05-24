 
getQT <- function(sig2, b, sig2.drift, b.drift, delta, driftMod)
{
	Qmat <- matrix(0, length(b), 3+2*driftMod)
	Tmat <- matrix(0, length(b), 2+2*driftMod)
	delta.sd <- max(getSD(delta), 16)
	if(!driftMod){
		idx <- (b<=0.1)
		psig2 <- sig2 #exp(log(sig2)-2*log(b))
		#

 	   Qmat[idx,1] <- sig2[idx]*ps1a(delta[idx],b[idx]) 
	   Qmat[!idx,1] <- psig2[!idx]*(delta[!idx] - 2*exp(pexp(delta[!idx],b[!idx],,TRUE)-log(b[!idx]))  + 
				 exp(pexp(delta[!idx],2*b[!idx],,TRUE)-log(2*b[!idx])))
						 
		#                             
		Qmat[idx,2] <- sig2[idx]*ps2a(delta[idx],b[idx])/2
		Qmat[!idx,2] <- psig2[!idx]*((1-2*exp(-b[!idx]*delta[!idx])+exp(-2*b[!idx]*delta[!idx]))/2)

		#
		Qmat[idx,3] <- sig2[idx]*ps3a(delta[idx],b[idx])
		#Qmat[!idx,3] <- sig2[!idx]*exp(pexp(delta[!idx],2*b[!idx],,TRUE)-log(2*b[!idx]))
		Qmat[!idx,3] <- sig2[!idx]*exp(log(b[!idx]) + pexp(delta[!idx],2*b[!idx],,TRUE))/2
		
		#
		Tmat[,1] <- exp(pexp(delta,b,,TRUE)-log(b))
		Tmat[,2] <- exp(-b*delta)
		
		#

# 		print(zapsmall(data.frame(v1=v1,v2=v2, Q=Qmat, xxx1, xxx2, xxx3, b)))
# 		cat(sum(Qmat[,1] - Qmat[,2]^2/Qmat[,3]<0), " bad values\n")
# 		cat(min(b),"\n")
# 		return(list(Qmat=Qmat, Tmat=Tmat))
	}
	else{
		psig2 <- sig2 #exp(log(sig2)-2*log(b))
		psig2.drift <- sig2.drift #exp(log(sig2.drift)-2*log(b.drift))
		V1 <- rep(0,length(b))
		V2 <- rep(0,length(b))
		idx <- (b<=0.01)
		idx2 <- (b.drift<=0.01)
		#
		V1[idx] <- sig2[idx]*ps1a(delta[idx],b[idx])
 		V1[!idx] <- psig2[!idx]*(delta[!idx] - 2*exp(pexp(delta[!idx],b[!idx],,TRUE)-log(b[!idx]))  + 
		                             exp(pexp(delta[!idx],2*b[!idx],,TRUE)-log(2*b[!idx])))
		#
		V2[idx2] <- sig2.drift[idx2]*ps1a(delta[idx2],b.drift[idx2])
		V2[!idx2] <- psig2.drift[!idx2]*(delta[!idx2] - 2*exp(pexp(delta[!idx2],b.drift[!idx2],,TRUE)-log(b.drift[!idx2]))  + 
		                             exp(pexp(delta[!idx2],2*b.drift[!idx2],,TRUE)-log(2*b.drift[!idx2])))
		#
	  	Qmat[,1] <- V1 + V2
	  	
	  	#
	  	Qmat[idx,2] <- sig2[idx]*ps2a(delta[idx],b[idx])/2
      	Qmat[!idx,2] <- psig2[!idx]*((1-2*exp(-b[!idx]*delta[!idx])+exp(-2*b[!idx]*delta[!idx]))/2)
      	#
	  	Qmat[idx2,3] <- sig2.drift[idx2]*ps2a(delta[idx2],b.drift[idx2])/2
      	Qmat[!idx2,3] <- psig2.drift[!idx2]*((1-2*exp(-b.drift[!idx2]*delta[!idx2])+exp(-2*b.drift[!idx2]*delta[!idx2]))/2)
      	#
		Qmat[idx,4] <- sig2[idx]*ps3a(delta[idx],b[idx])
		#Qmat[!idx,4] <- sig2[!idx]*exp(pexp(delta[!idx],2*b[!idx],,TRUE)-log(2*b[!idx]))
		Qmat[!idx,4] <- sig2[!idx]*exp(log(b[!idx]) + pexp(delta[!idx],2*b[!idx],,TRUE))/2
		#
		Qmat[idx2,5] <- sig2.drift[idx2]*ps3a(delta[idx2],b.drift[idx2])
		#Qmat[!idx2,5] <- sig2.drift[!idx2]*exp(pexp(delta[!idx2],2*b.drift[!idx2],,TRUE)-log(2*b.drift[!idx2]))
        Qmat[!idx2,5] <- sig2.drift[!idx2]*exp(log(b.drift[!idx]) + pexp(delta[!idx],2*b.drift[!idx],,TRUE))/2
		##
		##
		Tmat[idx,1] <- ps4(delta[idx],b[idx])
		Tmat[!idx,1] <- exp(pexp(delta[!idx],b[!idx],,TRUE)-log(b[!idx]))
		#
		Tmat[,2] <- exp(-b*delta)
		#
		Tmat[idx2,3] <- ps4(delta[idx2],b.drift[idx2])
		Tmat[!idx2,3] <- exp(pexp(delta[!idx2],b.drift[!idx2],,TRUE)-log(b.drift[!idx2]))
		#
		Tmat[,4] <- exp(-b.drift*delta)
	}
	Qmat <- round(Qmat,delta.sd)
	Tmat <- round(Tmat, delta.sd)
	return(list(Qmat=Qmat, Tmat=Tmat))
}

ps1 <- function(d,x){
	# Expansion of (d - 2*(1-exp(-d*x))/x + (1-exp(-2*d*x))/(2*x))/(x^2)
	d^3/3 - (d^4 * x)/4 + (7 * d^5 * x^2)/60 - (d^6 * x^3)/24 + (31 * d^7 * x^4)/2520 - 
	(d^8 * x^5)/320 + (127 * d^9 * x^6)/181440 - (17 * d^10 * x^7)/120960 + 
	(73 * d^11 * x^8)/2851200 - (31 * d^12 * x^9)/7257600 + (2047 * d^13 * x^10)/3113510400
}

ps1a <- function(d,x){
	# Expansion of d - 2*(1-exp(-d*x))/x + (1-exp(-2*d*x))/(2*x)
	(d^3 * x^2)/3 - (d^4 * x^3)/4 + (7 * d^5 * x^4)/60 - (d^6 * x^5)/24 + (31 * d^7 * x^6)/2520 - 
	(d^8 * x^7)/320 + (127 * d^9 * x^8)/181440 - (17 * d^10 * x^9)/120960 + (73 * d^11 * x^10)/2851200
}


ps2 <- function(d,x){
    # Expansion of (1 - 2*exp(-d*x) + exp(-2*d*x))/(x^2)
	d^2 - (d^3 * x) + (7 * d^4 * x^2)/12 - (d^5 * x^3)/4 + (31 * d^6 * x^4)/360 - 
	(d^7 * x^5)/40 + (127 * d^8 * x^6)/20160 - (17 * d^9 * x^7)/12096 + (73 * d^10 * x^8)/259200 - 
	(31 * d^11 * x^9)/604800 + (2047 * d^12 * x^10)/239500800
}

ps2a <- function(d,x){
	# Expansion of 1 - 2*exp(-d*x) + exp(-2*d*x)
	(d^2 * x^2) - (d^3 * x^3) + (7 * d^4 * x^4)/12 - (d^5 * x^5)/4 + (31 * d^6 * x^6)/360 - (d^7 * x^7)/40 + 
	(127 * d^8 * x^8)/20160 - (17 * d^9 * x^9)/12096 + (73 * d^10 * x^10)/259200
}

ps3 <- function(d,x){
    # Expansion of (1-exp(-2*d*x))/(2*x)   
	d - (d^2 * x) + (2 * d^3 * x^2)/3 - (d^4 * x^3)/3 + (2 * d^5 * x^4)/15 - (2 * d^6 * x^5)/45 +
	(4 * d^7 * x^6)/315 - (d^8 * x^7)/315 + (2 * d^9 * x^8)/2835 - (2 * d^10 * x^9)/14175 +
	(4 * d^11 * x^10)/155925 - (2 * d^12 * x^11)/467775 + (4 * d^13 * x^12)/6081075
}

ps3a <- function(d,x){
	# Expansion of x*(1-exp(-2*d*x))/2   
	(d * x^2) - (d^2 * x^3) + (2 * d^3 * x^4)/3 - (d^4 * x^5)/3 + (2 * d^5 * x^6)/15 - (2 * d^6 * x^7)/45 + 
	(4 * d^7 * x^8)/315 - (d^8 * x^9)/315 + (2 * d^9 * x^10)/2835
}

ps4 <- function(d,x){
	#Expansion of (1-exp(-d*b))/b
	d-(d^2 * x)/2+(d^3 * x^2)/6-(d^4 * x^3)/24+(d^5 * x^4)/120-
	(d^6 * x^5)/720+(d^7 * x^6)/5040-(d^8 * x^7)/40320+(d^9 * x^8)/362880-(d^10 * x^9)/3628800+
	(d^11 * x^10)/39916800-(d^12 * x^11)/479001600+(d^13 * x^12)/6227020800
}
