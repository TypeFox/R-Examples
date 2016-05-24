rLCD    <-	function(n, mle){
  
    phiknot		<-	mle$phi[mle$IsKnot==1]
    Fhatknot	<-	mle$Fhat[mle$IsKnot==1]
    xxknot		<-	mle$x[mle$IsKnot==1]
    alpha		<-	diff(phiknot)/diff(xxknot)
    beta		<-	phiknot[-length(phiknot)]-alpha*xxknot[-length(xxknot)]
    p		  	<-	diff(Fhatknot)
    beta2		<-	beta-log(p)
  
    k			<-	sample(1:length(p), n, replace=TRUE, prob=p)
    u			<-	runif(n)
    beta2k		<-	beta2[k]
    alphak		<-	alpha[k]
    xxk			<-	xxknot[k]
    res		  	<-	log(alphak*u*exp(-beta2k)+exp(alphak*xxk))/alphak
    return(res)
}