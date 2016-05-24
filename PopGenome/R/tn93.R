tn93 <- function(gT,gC,gG,gA,P1,P2,Q,length) {

#P1 is the A/G, P2 is T/C and Q is transversions
	
	sum <- gT+gC+gG+gA
	gT  <- gT/sum
	gC  <- gC/sum
	gG  <- gG/sum
	gA  <- gA/sum
	P1  <- P1/length
	P2  <- P2/length
	Q   <-  Q/length
	
	gR  <- gA + gG
	gY  <- gT + gC
	
        check <- (1-gR/(2.*gA*gG) * P1 - 1./(2.*gR) * Q)

        if(check<=0 | is.na(check) ){d1 <- NaN}else{
	   d1 <- -2.*gA*gG/gR * log(1. - gR/(2.*gA*gG) * P1 - 1./(2.*gR) * Q)
        }	

        check <- (1. - gY/(2.*gT*gC) * P2 - 1./(2.*gY) * Q) 
        if(check<=0 | is.na(check)){d2 <- NaN}else{
            d2 <- -2.*gT*gC/gY * log(1. - gY/(2.*gT*gC) * P2 - 1./(2.*gY) * Q)
        }
	
        check <- (1. - 1./(2.*gR*gY) * Q)
        if(check<=0 | is.na(check)){d3 <- NaN}else{
            d3 <- -2.*(gR*gY - gA*gG*gY/gR - gT*gC*gR/gY) * log(1. - 1./(2.*gR*gY) * Q)
	}

	d <- d1 + d2 + d3
      
      if(d=="NaN"){return(NaN)}
	if(d < 0){return (NaN)}
	return (d)
     
}
