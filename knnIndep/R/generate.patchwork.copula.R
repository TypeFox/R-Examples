generate.patchwork.copula <-
function(
		p = matrix(rbeta(bins*bins,alpha,beta),ncol=bins) #starting distribution of mass
		#parameters of beta distribution (if used)
		,alpha=.01
		,beta=1
		,c=1 #sharpness factors/concentration/dilution factor
		,npoints = 320,bins=20 #number of points to draw and grid size
		,returnmi=FALSE #return the mutual information along with the x and y data?
		,plot=FALSE #plot the data (does not open a new device)
){
	
	p = p^c/sum(p^c)
	
	nr = sample(0:(bins*bins-1),npoints,replace=TRUE,prob=p)
	i = (nr %% bins) 
	j = (nr %/% bins)
	
	pi = rowSums(p)
	pj = colSums(p)
	rescali=cumsum(pi)/sum(pi)
	rescalj=cumsum(pj)/sum(pj)
	mi = sum(p*log(t(t(p/pi)/pj)),na.rm=T)
	
	y = rescali[i+1] - runif(npoints)*pi[i+1]
	x = rescalj[j+1] - runif(npoints)*pj[j+1]
	
	if(plot){
		plot(x,y,main=mi)#,pch=".")
	}
	if(returnmi){
		return(list(x=x,y=y,mi=mi))
	}else{
		return(list(x=x,y=y))		
	}
	
}
