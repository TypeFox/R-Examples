#######################################################
#	ForwardSearch package
#	Bent Nielsen, 10 September 2014, version 1
#######################################################
#	Copyright 2014 Bent Nielsen
#	Nuffield College, OX1 1NF, UK
#	bent.nielsen@nuffield.ox.ac.uk
#
#	This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#######################################################

ForwardSearch.pointwise.asymptotics	<- function(psi,ref.dist="normal")
#	BN, 8 Sep 2014
#	This produces pointwise standard errors for Forward Search
#	Based on
#	Johansen & Nielsen (2013) Asymptotic analysis of the Forward Search
#	in:		psi				Number.  Takes value in interval 0,1
#			ref.dist		Character.  Reference distribution
#							"normal"	standard normal distribution
#	out:	varpi			Number.  sdv for X process
#			zeta			Number.	 consistency correction factor
#			sdv.unbiased	Number.  varpi/2/f
#			sdv.biased		Number.  varpi/2/f/zeta
#			c				Number.  c (median in unbiased case)
#			median.biased	Number.  median (in biased case)
{	#	ForwardSearch.pointwise.asymptotics
	if(ref.dist=="normal")
	{
		c			<- qnorm((psi+1)/2,0,1)
		f			<- dnorm(c,0,1)
		tau			<- psi - 2*c*f
		varkappa 	<- 3*psi - 2*c*(c^2+3)*f
	}
	wGG		<- psi*(1-psi)
	wGH 	<- tau*(1-psi)
	wHH		<- varkappa-tau^2
	a		<- 1 - c^3*f/tau
	b		<- c*f/tau
	w		<- a^2*wGG + 2*a*b*wGH + b^2*wHH
	sqrt.w	<- sqrt(w)	
	zeta	<- sqrt(tau/psi)
	sdv.unbiased	<- sqrt.w/2/f
	return(list(varpi			=sqrt.w				,
				zeta			=zeta				,
				sdv.unbiased	=sdv.unbiased		,
				sdv.biased		=sdv.unbiased/zeta	,
				c				=c					,
				median.biased	=c/zeta				))
}	#	ForwardSearch.pointwise.asymptotics

ForwardSearch.fit	<- function(x.1,y,psi.0=0.5,m.0=NULL,beta.0=NULL)
#	BN, 8 Sep 2014
#	This fits the Forward Search
#	Based on
#	Johansen & Nielsen (2013) Asymptotic analysis of the Forward Search
# 	dimensions:		n		Number of observations
#					dim.x	Number of regressors (including intercept)
#	in:		x.1				Matrix of dimension n x (dim.x -1).
#							Design matrix for regressors apart from constant. 
#			y				Vector of dimension n.
#							Dependent variable.
#			psi.0 			proportion of observations in initial set of set of selected observations.
#							Default is 0.5.
#							Initial set has round(n*psi.0) observations.
#			m.0				Number of observations in initial set of selected observations.
#							Default is NULL.
#							If value is given this overrides psi.0.
#			beta.0			Vector of dimension dim.x.
#							Initial estimator for regression coefficient.						
#							Default is NULL, which results in Least Trimmed Squares estimator
#							through beta.0	<- ltsReg(y~x.1,alpha=psi.0)$coefficients
#	out:	forward.beta			Matrix of dimension n x p.
#									Forward Search estimates of beta.		 
#			forward.sigma2.biased	Matrix of dimension	n x 1.
#									Forward Search estimates of sigma.
#									Values are *not* bias corrected.
#			forward.residual		Matrix of dimension	n x 1.
#									Forward Search estimates of forward residuals.
#									Values are *not* bias corrected.
#			m.0				Number of observations in initial set of selected observations.
#			y				Vector of dimension n.
#							Dependent variable from argument.
#			x				Matrix of dimension n x dim.x.
#							Design matrix for regressors.
#							Dependent variable from argument augmented with constant.
#							First column is constant.
{	#	ForwardSearch.fit
	########################
	#	get initial estimator
	if(is.null(beta.0)==TRUE)
		beta.0	<- ltsReg(y~x.1,alpha=psi.0)$coefficients
	########################
	#	get size of initial set
	if(is.null(psi.0)==TRUE & is.null(m.0)==TRUE )
	{
		psi.0	<- 0.5
		m.0		<- round(n*psi.0)
	}
	if(is.null(psi.0)==FALSE & is.null(m.0)==TRUE )
		m.0		<- round(n*psi.0)
	if(is.null(m.0)==FALSE )
		psi.0	<- m.0/n
	m.0	<- as.integer(m.0)	
	########################		
	#	get dimensions
	n		<- length(y)
	dim.x	<- length(beta.0)
	########################		
	#	get regressor matrix
	x	<- cbind(seq(1,1,length=n),x.1)
	########################		
	#	declare forward matrices
 	forward.beta			<- matrix(data=NA,n,dim.x)
	forward.residual		<- matrix(data=NA,n,1)
	forward.sigma2.biased	<- matrix(data=NA,n,1)
	forward.beta[m.0,]		<- beta.0
	########################		
	#	loop: forward search	
	for(m in m.0:(n-1))
	{
		abs.res	<- abs(y- x %*% forward.beta[m,] )
		ranks	<- sort(abs.res,index.return=TRUE)$ix
		forward.residual[m,1]	<-	abs.res[ranks[m+1],1]
		S		<- ranks[1:(m+1)]
		fit		<- glm.fit(x[S,],y[S])
		forward.beta[m+1,]		<- fit$coefficients
		forward.sigma2.biased[m+1,]	<- fit$deviance/(m+1)		#	glm.fit returns RSS as deviance in Gaussian case 
	}
	return(list(forward.beta			=forward.beta			,
				forward.sigma2.biased	=forward.sigma2.biased	,
				forward.residual		=forward.residual		,
				m.0						=m.0					,
				x						=x  	 				,
				y						=y    					
				))
}	#	ForwardSearch.fit

ForwardSearch.plot	<- function(FS,ref.dist="normal",bias.correct=FALSE,return=FALSE,plot.legend=TRUE,col=NULL,legend=NULL,lty=NULL,lwd=NULL,main=NULL,type=NULL,xlab=NULL,ylab=NULL)
#	BN, 8 Sep 2014
#	This gives the forward plot with simultaneous confidence bands
# 	based on Johansen and Nielsen (2013, 2014).
#	in:		FS				List. Value of the function ForwardSearch.fit
#			ref.dist		Character.  Reference distribution
#							"normal"	standard normal distribution
#			bias.correct	Logical.
#							If FALSE do not bias correct variance, so plots have
#							appearance similar to Atkinson and Riani (2000)
#							If TRUE do bias correct variance, so plots start at origin.
#							Default is FALSE
{	#	ForwardSearch.plot
	########################		
	#	get input
	m.0						<- FS$m.0
	dim.x					<- ncol(FS$x)
	n						<- nrow(FS$x)
	psi.0					<- m.0/n
	forward.residual		<- FS$forward.residual
	forward.sigma2.biased	<- FS$forward.sigma2.biased
	forward.residual.scaled	<- forward.residual / sqrt(forward.sigma2.biased)
	########################		
	#	get asymptotics
	psi			<- seq(1,n)/n
	FS.asymp	<- ForwardSearch.pointwise.asymptotics(psi,ref.dist)
	median		<- FS.asymp$median.biased
	sdv			<- FS.asymp$sdv.biased
	zeta		<- FS.asymp$zeta
	########################		
	#	bias correct
	if(bias.correct==TRUE)
	{
		forward.residual.scaled	<- forward.residual.scaled * zeta
		median					<- median * zeta
		sdv						<- sdv * zeta		
	}
	########################		
	#	get plot input
	if(is.null(col)==TRUE)
		col		<- c(6,5,4,3,2,1)
	if(is.null(legend)==TRUE)
		legend	<- c("gauge 0.001","gauge 0.005","gauge 0.01","gauge 0.05","gauge 0.10","pointwise median")
	if(is.null(lty)==TRUE)
		lty		<- c(1,1,1,1,1,1)
	if(is.null(lwd)==TRUE)
		lwd		<- c(1,1,1,1,1,1)
	if(is.null(main)==TRUE)
		main	<- bquote("Forward residual plot " ~ m[0] ~" =" ~.(m.0))
	if(is.null(type)==TRUE)
		type	<- "o"
	if(is.null(xlab)==TRUE)
		xlab	<- "m"
	if(is.null(ylab)==TRUE)
		ylab	<- "forward residual"
	########################		
	#	cut off values
	#	Taken from Table 3 of Johansen and Nielsen (2014)
	cut.off	<- matrix(nrow=6,ncol=10)
	cut.off[6,] <- c(0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   )
	cut.off[5,]	<- c(2.50,2.43,2.28,2.14,1.99,1.81,1.60,1.31,0.82,NA  )
	cut.off[4,]	<- c(2.77,2.71,2.58,2.46,2.33,2.19,2.02,1.79,1.45,0.69)
	cut.off[3,]	<- c(3.30,3.24,3.14,3.04,2.94,2.83,2.71,2.55,2.33,1.91)
	cut.off[2,]	<- c(3.49,3.44,3.35,3.26,3.15,3.04,2.95,2.81,2.62,2.26)
	cut.off[1,]	<- c(3.90,3.85,3.77,3.69,3.62,3.53,3.43,3.32,3.18,2.92)
	colnames(cut.off) <- c("0.05","0.10","0.20","0.30","0.40","0.50","0.60","0.70","0.80","0.90")
	rownames(cut.off) <- legend
	########################		
	#	plot
	m	<- seq(from=m.0,to=n-1,by=1)
	plot(m,forward.residual.scaled[m],type=type,xlab=xlab,ylab=ylab,main=main)
	n.cut.off	<- nrow(cut.off)
	if(psi.0<=1.00 && psi.0>0.85)	cut.off.choice	<- 10
	if(psi.0<=0.85 && psi.0>0.75)	cut.off.choice	<- 9
	if(psi.0<=0.75 && psi.0>0.65)	cut.off.choice	<- 8
	if(psi.0<=0.65 && psi.0>0.55)	cut.off.choice	<- 7
	if(psi.0<=0.55 && psi.0>0.45)	cut.off.choice	<- 6
	if(psi.0<=0.45 && psi.0>0.35)	cut.off.choice	<- 5
	if(psi.0<=0.35 && psi.0>0.25)	cut.off.choice	<- 4
	if(psi.0<=0.25 && psi.0>0.15)	cut.off.choice	<- 3
	if(psi.0<=0.15 && psi.0>0.075)	cut.off.choice	<- 2
	if(psi.0<=0.075)  				cut.off.choice	<- 1
	for(i in 1:n.cut.off)
		lines(m,median[m]+cut.off[i,cut.off.choice]*sdv[m]/sqrt(n-dim.x),lty=lty[i],col=col[i],lwd=lwd[i])
	if(plot.legend==TRUE)	
		legend(x="topleft",legend=legend,lty=lty,col=col,lwd=lwd)
	########################		
	#	return
	if(return==TRUE)
		return(list(ref.dist				= ref.dist							,
					bias.correct			= bias.correct						,
					forward.residual.scaled	= forward.residual.scaled * zeta	,
					forward.asymp.median	= median * zeta						,
					forward.asymp.sdv		= sdv * zeta						,
					cut.off					= cut.off
					))		
}	#	ForwardSearch.plot

ForwardSearch.stopped	<- function(FS,m)
#	BN, 8 Sep 2014
#	A Forward Search gives a sequence of regression estimators
#	This function gives the regression estimators when stopped at m
#	in:		FS				List. Value of the function ForwardSearch.fit
#			m				Integer. Stopping time
#	out:	ranks.selected	Vector. Ranks of m observations in the selected set. 
#			ranks.outliers	Vector. Ranks of n-m observations that are not selected.
#							These are the "outliers". It is the complement to ranks.selected.
#			beta.m			Vector. Least squares estimator based on ranks.selected.
#			sigma2.biased	Scalar. Least squares residual variance based on ranks.selected.
#							Value is *not* bias corrected.
{	#	ForwardSearch.stopped
	########################		
	#	check m
	if(m <= FS$m.0)
		print("WARNING in ForwardSearch.stopped: m <= m.0")
	########################		
	#	get input
	n	<- length(FS$y)
	y	<- FS$y
	x	<- FS$x
	beta.previous	<- FS$forward.beta[m-1,]
	########################		
	#	compute residuals and their ranks
	res		<- y- x %*% beta.previous 
	abs.res	<- abs(res )
	ranks	<- sort(abs.res,index.return=TRUE)$ix
	########################
	#	return	
	return(list(res				=res										,
				ranks.selected	=ranks[1:m]									,
				ranks.outliers	=ranks[(m+1):n]								,
				beta.m			=FS$forward.beta[m,]						,
				sigma2.biased	=FS$forward.sigma2.biased[m,]				
				))
}	#	ForwardSearch.stopped


