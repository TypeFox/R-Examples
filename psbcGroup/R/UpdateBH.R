UpdateBH <-
function(survObj, priorPara, ini){
	
	n <- survObj$n
	p <- survObj$p	
	
	xbeta		<- ini$xbeta
	ind.r_d		<- priorPara$ind.r_d
	c0			<- priorPara$c0
	hPriorSh	<- priorPara$hPriorSh	
	d			<- priorPara$d	
	J			<- priorPara$J
		
	exp.xbeta	<- exp(xbeta)
	exp.xbeta.mat	<- matrix(rep(exp.xbeta, J), n, J)
		
	h.rate	<-  colSums(exp.xbeta.mat * ind.r_d) + c0
	h		<- rgamma(J, shape = hPriorSh + d, rate = h.rate)
	
	return(h)
	
	}

