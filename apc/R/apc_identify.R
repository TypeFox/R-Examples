#######################################################
#	apc package
#	Bent Nielsen, 3 April 2015, version 1.1
#	functions to identify parameters
#######################################################
#	Copyright 2014,2015 Bent Nielsen
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

#########################################################
#	apc.identify
#########################################################
apc.identify	<- function(apc.fit.model)
#	BN  2 feb 2016	Changed: parameter label: date to character changed to allow nice decimal points
#					using apc.internal.function.date.2.character
#	3 Mar 2015
#	In:		apc.fit.model			list		from apc.fit.model
#	Out:	index.age.max	 		vector. 	Indices for identified age parameters 
#			index.per.max			vector. 	Indices for identified period parameters 
#			index.coh.max			vector. 	Indices for identified cohort parameters 
#			dates.max				vector.		Dates for parameters
#			coefficients.ssdd		matrix.		4 columns: est, sdv, t, p
#												double sums of double differences
#			covariance.ssdd			matrix.		covariance matrix
#			coefficients.detrend	matrix.		4 columns: est, sdv, t, p         
#												detrended double sums of double differences 
#			covariance.detrend		matrix.		covariance matrix                 
{	#	function.apc.identify.double.sums
	##############################
	#	create indicator for mixed parametrisation
	mixed.par	<- isTRUE(apc.fit.model$model.family %in% c("poisson.response"))
	##############################
	#	model design list
	model.design.list.sub	<- c("AP","AC","PC") #,"Ad","Pd","Cd","A","P","C","t","tA","tP","tC","1")
	##############################
	#	get values
	coefficients	<- apc.fit.model$coefficients.canonical	#	4 columns
	covariance		<- apc.fit.model$covariance.canonical
	slopes			<- apc.fit.model$slopes
	difdif			<- apc.fit.model$difdif
	dates			<- apc.fit.model$dates
	index.age		<- apc.fit.model$index.age
	index.per		<- apc.fit.model$index.per
	index.coh		<- apc.fit.model$index.coh
	age1			<- apc.fit.model$age1
	per1			<- apc.fit.model$per1
	coh1			<- apc.fit.model$coh1
	unit			<- apc.fit.model$unit
	per.zero		<- apc.fit.model$per.zero
	per.odd			<- apc.fit.model$per.odd
	U				<- apc.fit.model$U
	age.max			<- apc.fit.model$age.max
	per.max			<- apc.fit.model$per.max
	coh.max			<- apc.fit.model$coh.max
	model.design	<- apc.fit.model$model.design
	n.decimal		<- apc.fit.model$n.decimal
	##############################
	#	derived values
	det.max	<- 1+sum(slopes)
	det.sub	<- 1+sum(slopes)-sum(difdif)
	xi.max	<- det.max+difdif[1]*age.max+difdif[2]*per.max+difdif[3]*coh.max
	xi.sub	<- det.sub+difdif[1]*age.max+difdif[2]*per.max+difdif[3]*coh.max
	xi.dif	<- det.sub+difdif[1]*(age.max-1)+difdif[2]*(per.max-1)+difdif[3]*(coh.max-1)
	xi		<- det.max+difdif[1]*(age.max-2)+difdif[2]*(per.max-2)+difdif[3]*(coh.max-2)
	##############################
	#	construct for indices for double sums of double difference parameters
	#	first for use with (double difference) canonical parameters
	index.age.max	<- NULL
	index.per.max	<- NULL
	index.coh.max	<- NULL
	start		<- det.max
	if(difdif[1])	{	index.age.max	<- start+seq(1,age.max);	start	<- start+age.max	}
	if(difdif[2])	{	index.per.max	<- start+seq(1,per.max);	start	<- start+per.max	}
	if(difdif[3])	{	index.coh.max	<- start+seq(1,coh.max);	start	<- start+coh.max	}
	#	then for use with in reparametrised submodels in terms of sums of differences 
	index.age.sub	<- NULL
	index.per.sub	<- NULL
	index.coh.sub	<- NULL
	if(model.design != "APC")
	{
		start		<- det.sub
		if(difdif[1])	{	index.age.sub	<- start+seq(1,age.max);	start	<- start+age.max	}
		if(difdif[2])	{	index.per.sub	<- start+seq(1,per.max);	start	<- start+per.max	}
		if(difdif[3])	{	index.coh.sub	<- start+seq(1,coh.max);	start	<- start+coh.max	}
	}
	#	then for use with in reparametrised submodels in terms of differences 
	index.age.dif	<- NULL
	index.per.dif	<- NULL
	index.coh.dif	<- NULL
	if(model.design != "APC")
	{
		start		<- det.sub
		if(difdif[1])	{	index.age.dif	<- start+seq(1,age.max-1);	start	<- start+age.max-1	}
		if(difdif[2])	{	index.per.dif	<- start+seq(1,per.max-1);	start	<- start+per.max-1	}
		if(difdif[3])	{	index.coh.dif	<- start+seq(1,coh.max-1);	start	<- start+coh.max-1	}
	}
	##############################
	#	construct dates for for double difference parameters 
	#	first for use with (double difference) canonical parameters
	dates.max		<- matrix(data=NA,nrow=xi.max,ncol=1)
	if(difdif[1])	dates.max[index.age.max,1]	<- age1+seq(0,age.max-1)*unit	
	if(difdif[2])	dates.max[index.per.max,1]	<- per1+seq(0,per.max-1)*unit
	if(difdif[3])	dates.max[index.coh.max,1]	<- coh1+seq(0,coh.max-1)*unit
	#	then for use with in reparametrised submodels of sums of differences 
	dates.sub		<- NULL
	if(model.design != "APC")
	{
		dates.sub		<- matrix(data=NA,nrow=xi.sub,ncol=1)
		if(difdif[1])	dates.sub[index.age.sub,1]	<- age1+seq(0,age.max-1)*unit	
		if(difdif[2])	dates.sub[index.per.sub,1]	<- per1+seq(0,per.max-1)*unit
		if(difdif[3])	dates.sub[index.coh.sub,1]	<- coh1+seq(0,coh.max-1)*unit
	}	
	#	then for use with in reparametrised submodels of sums of differences 
	dates.dif		<- NULL
	if(model.design != "APC")
	{
		dates.dif		<- matrix(data=NA,nrow=xi.dif,ncol=1)
		if(difdif[1])	dates.dif[index.age.dif,1]	<- age1+seq(1,age.max-1)*unit	
		if(difdif[2])	dates.dif[index.per.dif,1]	<- per1+seq(1,per.max-1)*unit
		if(difdif[3])	dates.dif[index.coh.dif,1]	<- coh1+seq(1,coh.max-1)*unit
	}	
	##############################
	#	get linear transformation matrix
	#		from canonical parameter
	#		to standard representation
	#	level + slope_age (i-1) + slope_coh (k-1)
	#		+ sum sum DD age 	[padded with zeros]
	#		+ sum sum DD period	[padded with zeros]
	#		+ sum sum DD cohort	[padded with zeros]
	#	a summation function is needed
	#	for sum sum DD age:		use with U=U
	#	for sum sum DD cohort:	use with U=U
	#	for sum sum DD period, L=per.odd=TRUE:	use with U=2
	#	for sum sum DD period, L=per.odd=FALSE:	use with U=1	
	function.ssdd	<- function(n,U)
	#	BN, 4 mar 2015
	#	U is the anchoring point in the summation
	{	#	function.ssdd
		m	<- matrix(data=0,nrow=n+2,ncol=n)
		if(U>1)
			for(row in 1:(U-1))
				m[row,row:(U-1)]	<- 1:(U-row)	
		if(U<n+1)		
			for(row in (U+2):(n+2))
				m[row,U:(row-2)]	<- (row-U-1):1	
		return(m)	
	}	#	function.ssdd
	##############################
	#	declare linear transformation matrix
	m.ssdd	<-	matrix(data=0,nrow=xi.max,ncol=xi)									
	m.ssdd[1:det.max,1:det.max]	<- diag(det.max)							#	level/trend terms
	if(difdif[1])	m.ssdd[index.age.max,index.age]	<- function.ssdd(age.max-2,U)			#	alpha
	if(difdif[2])	m.ssdd[index.per.max,index.per]	<- function.ssdd(per.max-2,per.odd+1)	#	beta 
	if(difdif[3])	m.ssdd[index.coh.max,index.coh]	<- function.ssdd(coh.max-2,U)			# 	gamma
	##############################
	#	get linear transformation matrix
	#		from standard representation
	#		to detrended representation
	#	this matrix more complicated because
	#	linear trends are moved from the
	#	time effects to the slopes.
	#
	#	a detrending function is needed
	function.detrend	<- function(n)
	#	BN, 3 apr 2015
	#	in:		n			is the dimension
	#	Out		m			matrix of dimension n x n
	#						for detrending an n vector,
	#						takes an identity matrix
	#						replaces first column with (col-n)/(n-1)
	#						replaces  last column with (1-col)/(n-1)
	{	#	function.detrend
		#	m defines the detrending
		m			<- diag(n);
		m[1:n,1]	<- (seq(1:n)-n)/(n-1);
		m[1:n,n]	<- (1-seq(1:n))/(n-1);
		m[1,1]		<- 0;
		m[n,n]		<- 0;
		return(m)
	}	#	function.detrend
	#	declare linear transformation matrix
	m.detrend	<-	diag(xi.max)
	#	move anchoring of linear trend from U to 1.
	if(sum(slopes)==2)
		m.detrend[1,2:3]	<- 1-U
	if(sum(slopes)==1 && !slopes[2])
		m.detrend[1,2]		<- 1-U
	#	detrend age effects, move linear trend to deterministics
	if(difdif[1])
	{	m.detrend[1,index.age.max[1]]	<- 1
		m.detrend[2,index.age.max[c(1,age.max)]]	<- c(-1,1)/(age.max-1)
		m.detrend[index.age.max,index.age.max]	<- function.detrend(age.max)		
	}
	if(difdif[3])
	{	# there are 2 slopes if slopes=c(1,0,1)
		# there is  1 slope  if slopes=c(0,0,1)
		# recall det.max <- 1+sum(slopes)
		m.detrend[1,index.coh.max[1]]	<- 1
		m.detrend[det.max,index.coh.max[c(1,coh.max)]]	<- c(-1,1)/(coh.max-1)
		m.detrend[index.coh.max,index.coh.max]	<- function.detrend(coh.max)		
	}
	if(difdif[2])
	{	# if slopes=c(1,0,1) the period slope gives age & cohort slopes with equal weight
		# if slopes=c(0,1,0) the period slope gives a period     slope
		if(!slopes[2])  m.detrend[1,index.per.max[c(1,per.max)]] <- c(1,0)+c(1,-1)*per.zero/(per.max-1)
		if(slopes[2]) 	m.detrend[1,index.per.max[c(1,per.max)]] <- c(1,0)
		if(slopes[1])	m.detrend[2,index.per.max[c(1,per.max)]]  <- c(-1,1)/(per.max-1)
		if(slopes[2])	m.detrend[2,index.per.max[c(1,per.max)]]  <- c(-1,1)/(per.max-1)
		if(slopes[3])	m.detrend[3,index.per.max[c(1,per.max)]]  <- c(-1,1)/(per.max-1)
		m.detrend[index.per.max,index.per.max]	<- function.detrend(per.max)		
	}
	##############################	
	##############################
	#	Now, manipulate estimates using m.ssdd and m.detrend
	##############################
	#	get estimates
	coefficients.ssdd		<- m.ssdd 		%*% coefficients
	coefficients.detrend	<- m.detrend 	%*% coefficients.ssdd
	##############################
	#	construct row names 
	names.ssdd		<- c("level")
	names.detrend	<- c("level")
	if(slopes[1])	{	names.ssdd		<- c(names.ssdd	  ,"age slope")
						names.detrend	<- c(names.detrend,"age slope")		}		
	if(slopes[2])	{	names.ssdd		<- c(names.ssdd	  ,"period slope")
						names.detrend	<- c(names.detrend,"period slope")	}
	if(slopes[3])	{	names.ssdd		<- c(names.ssdd	  ,"cohort slope")
						names.detrend	<- c(names.detrend,"cohort slope")	}		
	if(difdif[1])
		for(i in 1:age.max)
		{	names.ssdd		<- c(names.ssdd	  ,paste("SS_DD_age_"   ,apc.internal.function.date.2.character((dates.max[index.age.max,1])[i],n.decimal),sep=""))
			names.detrend	<- c(names.detrend,paste("SS_DD_age_"   ,apc.internal.function.date.2.character((dates.max[index.age.max,1])[i],n.decimal),sep=""))	}
	if(difdif[2])															      
		for(i in 1:per.max)															      
		{	names.ssdd		<- c(names.ssdd	  ,paste("SS_DD_period_",apc.internal.function.date.2.character((dates.max[index.per.max,1])[i],n.decimal),sep=""))
			names.detrend	<- c(names.detrend,paste("SS_DD_period_",apc.internal.function.date.2.character((dates.max[index.per.max,1])[i],n.decimal),sep=""))	}		
	if(difdif[3])																	      
		for(i in 1:coh.max)															      
		{	names.ssdd		<- c(names.ssdd	  ,paste("SS_DD_cohort_",apc.internal.function.date.2.character((dates.max[index.coh.max,1])[i],n.decimal),sep=""))
			names.detrend	<- c(names.detrend,paste("SS_DD_cohort_",apc.internal.function.date.2.character((dates.max[index.coh.max,1])[i],n.decimal),sep="")) }
	rownames(coefficients.ssdd	 )	<- names.ssdd		
	rownames(coefficients.detrend)	<- names.detrend
	##############################
	#	get covariance matrix, noting that if using a mixed parametrisation
	#	top right block of m should be zero, to reflect that intercept is not changed
	if(mixed.par)
	{	m.ssdd[1,2:xi]			<- 0
		m.detrend[1,2:xi.max]	<- 0
		}
	covariance.ssdd		<- m.ssdd 	%*% covariance 			%*% t(m.ssdd   )
	covariance.detrend	<- m.detrend%*% covariance.ssdd 	%*% t(m.detrend)
	##############################
	#	get standard errors 
	coefficients.ssdd[,2]		<- sqrt(diag(covariance.ssdd   ))
	coefficients.detrend[,2]	<- sqrt(diag(covariance.detrend))
	#	set NA for ad hoc identified entries
	if(difdif[1])	coefficients.ssdd[   index.age.max,2][U:(U+1)]	<- NA
	if(difdif[2])	coefficients.ssdd[   index.per.max,2][(per.odd+1):(per.odd+2)]	<- NA
	if(difdif[3])	coefficients.ssdd[   index.coh.max,2][U:(U+1)]	<- NA
	if(difdif[1])	coefficients.detrend[index.age.max,2][c(1,age.max)]	<- NA                           
	if(difdif[2])	coefficients.detrend[index.per.max,2][c(1,per.max)]	<- NA   
	if(difdif[3])	coefficients.detrend[index.coh.max,2][c(1,coh.max)]	<- NA                           
	######################
	#	get t-statistics & p-values
	function.get.t_and_p	<- function(coefficients)
	{
		coefficients[,3]	<- coefficients[,1] / coefficients[,2]
		coefficients[,4]	<- 2*pnorm(abs(coefficients[,3]),lower.tail=FALSE)
		return(coefficients)
	}
	coefficients.ssdd		<- function.get.t_and_p(coefficients.ssdd	)
	coefficients.detrend	<- function.get.t_and_p(coefficients.detrend)
	#################################
	#################################
	#	Now turn to demean and dif
	##############################
	#	get linear transformation matrix
	#		from detrended representation
	#		to demeaned levels
	#############################
	#	default values if model design not right
	coefficients.demean	<- NULL
	covariance.demean	<- NULL
	coefficients.dif	<- NULL
	covariance.dif		<- NULL
	if(isTRUE(model.design %in% model.design.list.sub)==TRUE)
	{	################################
		#	linear transformation for demean
		m.demean		<- matrix(data=0,nrow=xi.sub,ncol=xi.max)
		m.demean[1,1]	<- 1
		m.demean[(det.sub+1):xi.sub,(det.max+1):xi.max]	<- diag(xi.sub-det.sub)
		if(model.design=="AC")
		{	m.demean[index.age.sub,2]	<- seq(0,age.max-1)
			m.demean[index.coh.sub,3]	<- seq(0,coh.max-1)
		}
		if(model.design=="AP")
		{	m.demean[index.age.sub,2]	<- seq(0,age.max-1)
			m.demean[index.age.sub,3]	<- seq(0,age.max-1)*(-1)
			m.demean[index.per.sub,3]	<- seq(0,per.max-1)
			m.demean[1,3]				<- age.max-1
		}
		if(model.design=="PC")
		{	m.demean[index.per.sub,2]	<- seq(0,per.max-1)
			m.demean[index.coh.sub,2]	<- seq(0,coh.max-1)*(-1)
			m.demean[1,2]				<- age.max-1
			m.demean[index.coh.sub,3]	<- seq(0,coh.max-1)
		}
		##############################
		#	get linear transformation matrix for dif
		function.dif	<- function(n)
		##############################
		#	BN, 3 dec 2013
		#	in:		n			is the dimension of the block
		#	Out		m			matrix dimension n-1 x n
		#						for difference an n-vector
		{	#	function.dif
			m	<- diag(n)
			m[2:n,1:(n-1)]	<- m[2:n,1:(n-1)] - diag(n-1)	
			return(m[2:n,1:n])
		}	#	function.dif
		#	declare transformation
		m.dif		<- NULL
		m.dif			<- matrix(data=0,nrow=xi.dif,ncol=xi.sub)
		m.dif[1:det.sub,1:det.sub]			<- diag(det.sub)
		m.dif[index.age.dif,index.age.sub]	<- function.dif(age.max)
		m.dif[index.per.dif,index.per.sub]	<- function.dif(per.max)
		m.dif[index.coh.dif,index.coh.sub]	<- function.dif(coh.max)
		##############################	
		#	Now, manipulate estimates using m.demean and m.dif
		##############################
		#	get estimates
		coefficients.demean		<- m.demean 	%*% coefficients.detrend
		coefficients.dif		<- m.dif 		%*% coefficients.demean
		##############################
		#	construct row names 
		names.demean	<- c("level")
		names.dif		<- c("level")
		if(difdif[1])
			for(i in 1:age.max)
			{	names.demean	<- c(names.demean ,paste(  "S_D_age_"   ,apc.internal.function.date.2.character((dates.sub[index.age.sub,1])[i],n.decimal),sep=""))
				if(i>1)
				names.dif		<- c(names.dif	  ,paste(    "D_age_"   ,apc.internal.function.date.2.character((dates.sub[index.age.sub,1])[i],n.decimal),sep=""))	}
		if(difdif[2])															      
			for(i in 1:per.max)															      
			{	names.demean	<- c(names.demean ,paste(  "S_D_period_",apc.internal.function.date.2.character((dates.sub[index.per.sub,1])[i],n.decimal),sep=""))
				if(i>1)
				names.dif		<- c(names.dif	  ,paste(    "D_period_",apc.internal.function.date.2.character((dates.sub[index.per.sub,1])[i],n.decimal),sep=""))	}		
		if(difdif[3])																	      
			for(i in 1:coh.max)															      
			{	names.demean	<- c(names.demean ,paste(  "S_D_cohort_",apc.internal.function.date.2.character((dates.sub[index.coh.sub,1])[i],n.decimal),sep=""))
				if(i>1)
				names.dif		<- c(names.dif	  ,paste(    "D_cohort_",apc.internal.function.date.2.character((dates.sub[index.coh.sub,1])[i],n.decimal),sep=""))	}
		rownames(coefficients.demean )	<- names.demean	
		rownames(coefficients.dif	 )	<- names.dif	
		##############################
		#	get covariance matrix, noting that if using a mixed parametrisation
		#	top right block of m should be zero, to reflect that intercept is not changed
		if(mixed.par)
		{	m.demean[1,2:xi.sub]	<- 0
			m.dif[1,2:xi.sub]		<- 0
			}
		covariance.demean	<- m.demean %*% covariance.detrend	%*% t(m.demean )
		covariance.dif		<- m.dif	%*% covariance.demean	%*% t(m.dif	   )
		##############################
		#	get standard errors 
		coefficients.demean[,2]		<- sqrt(diag(covariance.demean ))
		coefficients.dif[,2]		<- sqrt(diag(covariance.dif	   ))
		#	set NA for ad hoc identified entries
		if(difdif[1])	coefficients.demean[ index.age.sub,2][1]	<- NA                           
		if(difdif[2])	coefficients.demean[ index.per.sub,2][1]	<- NA   
		if(difdif[3])	coefficients.demean[ index.coh.sub,2][1]	<- NA
		coefficients.demean 	<- function.get.t_and_p(coefficients.demean )
		coefficients.dif    	<- function.get.t_and_p(coefficients.dif    )
	}	
	##############################
	return(list(index.age.max	 		= index.age.max    		,
				index.per.max			= index.per.max    		,
				index.coh.max			= index.coh.max    		,
				dates.max				= dates.max        		,
				index.age.sub	 		= index.age.sub	  		, 
				index.per.sub			= index.per.sub	  		, 
				index.coh.sub			= index.coh.sub	  		, 
				dates.sub				= dates.sub		  		, 
				index.age.dif	 		= index.age.dif	 		, 
				index.per.dif			= index.per.dif	 		, 
				index.coh.dif			= index.coh.dif	 		, 
				dates.dif				= dates.dif		 		, 
				coefficients.ssdd		= coefficients.ssdd		,
				covariance.ssdd		 	= covariance.ssdd	 	,
				coefficients.detrend	= coefficients.detrend	,
				covariance.detrend		= covariance.detrend	, 
				coefficients.demean		= coefficients.demean	,
				covariance.demean	 	= covariance.demean		,
				coefficients.dif		= coefficients.dif		,
				covariance.dif			= covariance.dif		 
		))
}	#	apc.identify

