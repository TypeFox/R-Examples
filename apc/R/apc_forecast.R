#######################################################
#	apc package
#	Bent Nielsen, 8 Jan 2016, version 1.2
#	Forecasting
#######################################################
#	Copyright 2016 Bent Nielsen
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

apc.forecast.ac <- function(apc.fit,sum.per.by.age=NULL,sum.per.by.coh=NULL,covariance.output=FALSE,suppress.warning=TRUE)
#	BN 31 Jan 2016
#	Based on Martinez Miranda, Nielsen and Nielsen (2015) JRSS-A
#	Requires AC structure
#	Computes linear predictors for all models
#	Computes distribution forecasts for a Poisson response model.
#	This is done for the triangle which shares age and cohort indices with the data.
#	in
#			apc.fit						List. Output for apc.fit.model
#										Note: apc.fit.model should be run first for a Poisson response model with AC structure so that
#										apc.fit$model.design=="AC"
#										apc.fit$model.family=="poisson.response"
#			sum.per.by.age
#			sum.per.by.coh
#			covariance.output			Logical. If TRUE gives output for computing the covariance matrix for
#										cell-by-cell point forecasts. This matrix can be very large.
#	out
#			linear.predictors.forecast	Vector. Linear predictors for forecast area
#			index.trap.J				Matrix. age-coh coordinates for vector. Similar structure to index.trap in apc.index
#			trap.response.forecast		Matrix. Includes data and point forecasts. Forecasts in lower right triangle.
#			response.forecast.cell		Matrix. 4 columns.
#										1: Point forecasts.
#										2: corresponding forecast standard errors
#										3: process standard errors
#										4: estimation standard errors
#										Note that the square of column 2 equals the sums of squares of columns 3 and 4
#										Note that index.trap.J gives the age-coh coordinates for each entry 
#			response.forecast.age		Same as response.forecast.cell, but point forecasts cumulated by age.
#			response.forecast.per		Same as response.forecast.cell, but point forecasts cumulated by per.
#			response.forecast.coh		Same as response.forecast.cell, but point forecasts cumulated by coh.
#			response.forecast.all		Same as response.forecast.cell, but point forecasts cumulated by age & coh.
#			response.forecast.per.ic	Same as response.forecast.cell, but point forecasts cumulated by per and intercept corrected by
#										multiplying column 1 of response.forecast.per by intercept.correction.per
#			intercept.correction.per	Numeric.
#										The intercept correction is constructed as the ratio of
#										the sum of data entries for the last period and the sum of the corresponding fitted values
#			covariance.proc				Matrix. Only if covariance.output=TRUE.
#			covariance.est				Vector. Only if covariance.output=TRUE. 
#										Note overall covariance matrix given by
#										covariance.forecast=covariance.proc+diag(covariance.est)
{	#	apc.forecast.ac
	########################
	#	get values
	age1			        <- apc.fit$age1			   	
	per1			        <- apc.fit$per1
	coh1			        <- apc.fit$coh1
	unit					<- apc.fit$unit
	age.max			        <- apc.fit$age.max			   	
	per.max			        <- apc.fit$per.max	
	coh.max			        <- apc.fit$coh.max	
	per.odd			        <- apc.fit$per.odd
	per.zero		        <- apc.fit$per.zero
	U				        <- apc.fit$U
	n.data					<- apc.fit$n.data
	response		        <- apc.fit$response
	index.data		        <- apc.fit$index.data
	index.trap		        <- apc.fit$index.trap
	model.design			<- apc.fit$model.design
	model.family			<- apc.fit$model.family
	coefficients.canonical	<- apc.fit$coefficients.canonical
	covariance.canonical	<- apc.fit$covariance.canonical
	fitted.values			<- apc.fit$fitted.values
	data.format				<- apc.fit$data.format
	n.decimal				<- apc.fit$n.decimal
	#	derived values
	per.forecast.J			<- age.max+coh.max-1-per.zero-per.max
	n.data.J				<- (per.forecast.J*(per.forecast.J+1)) %/% 2	# note %/% has higher precedence than *
	##############################
	#	Check Input
	if(isTRUE(model.design=="AC")==FALSE)
		return(cat("ERROR apc.forecast.ac: model.design must be 'AC'\n"))
	if(per.forecast.J==0)
		return(cat("ERROR apc.forecast.ac: forecast area empty since per.forecast=age.max+coh.max-1-per.zero-per.max=0\n"))
	distribution.forecasts	<- TRUE	
	if(isTRUE(model.family=="poisson.response")==FALSE)
	{	return(cat("WARNING apc.forecast.ac: model.family is not 'poisson.response': only point forecasts are generated\n"))
		distribution.forecasts	<- FALSE
	}
	if(is.null(sum.per.by.age)==FALSE && distribution.forecasts==TRUE)
	{	if(is.vector(sum.per.by.age)==FALSE) 
			return(cat("ERROR apc.forecast.ac: sum.per.by.age argument should be vector\n"))
		if(length(sum.per.by.age)==1)	
			sum.per.by.age <- c(1,1)*sum.per.by.age
		if(sum.per.by.age[1]>sum.per.by.age[2])
			return(cat("ERROR apc.forecast.ac: sum.per.by.age argument should have element 1 <= element 2'\n"))
		if(sum.per.by.age[1]<age.max-per.forecast.J+1)
			return(cat("ERROR apc.forecast.ac: sum.per.by.age[1] too small'\n"))
		if(sum.per.by.age[2]>age.max)	
			return(cat("ERROR apc.forecast.ac: sum.per.by.age[2] too large'\n"))
	}		
	if(is.null(sum.per.by.coh)==FALSE && distribution.forecasts==TRUE)
	{	if(is.vector(sum.per.by.coh)==FALSE) 
			return(cat("ERROR apc.forecast.ac: sum.per.by.coh argument should be vector\n"))
		if(length(sum.per.by.coh)==1)	
			sum.per.by.coh <- c(1,1)*sum.per.by.coh
		if(sum.per.by.coh[1]>sum.per.by.coh[2])
			return(cat("ERROR apc.forecast.ac: sum.per.by.coh argument should have element 1 <= element 2'\n"))
		if(sum.per.by.coh[1]<coh.max-per.forecast.J+1)
			return(cat("ERROR apc.forecast.ac: sum.per.by.coh[1] too small'\n"))
		if(sum.per.by.coh[2]>coh.max)	
			return(cat("ERROR apc.forecast.ac: sum.per.by.coh[2] too large'\n"))
	}		
	##############################
	#	Generate forecast trapezoid
	trap.response.forecast	<- apc.data.list.subset(apc.fit,suppress.warning=suppress.warning)$response
	#	index for forecast area
	index.trap.J		<- matrix(nrow=n.data.J,ncol=2)
	s	<- 0
	for(row in 1:per.forecast.J)
	{	index.trap.J[(s+1):(s+row),1]	<- seq(age.max-row+1,age.max)	# age
		index.trap.J[(s+1):(s+row),2]	<- coh.max-per.forecast.J+row	# cohort
		s	<- s + row
	}
	colnames(index.trap.J)	<- c("age","coh")
	##############################
	#	Get design for forecast area
	#	apc index list for forecast area
	apc.index.J	<- list(	age.max		=  age.max			,		
							per.max		=  age.max+coh.max-1,
							coh.max		=  coh.max			,
							index.trap	=  index.trap.J		,
							per.zero	=  per.zero			,
							per.odd		=  per.odd			,
							U			=  U				,
							n.data		=  n.data.J			)
	#	get design matrix
   	design.J	<- apc.get.design(apc.index.J,model.design)$design
	##############################
	#	Get design for sample area
	#	apc index list for sample area
	apc.index	<- list(	age.max		=  age.max			,		
							per.max		=  per.max			,
							coh.max		=  coh.max			,
							age1		=  age1				,	# this entry used in intercept correction
							per1		=  per1				,	# this entry used in intercept correction
							coh1		=  coh1				,	# this entry used in intercept correction
							unit		=  unit				,	# this entry used in intercept correction
							index.trap	=  index.trap 		,
							index.data	=  index.data		,	# this entry used in intercept correction
							per.zero	=  per.zero			,
							per.odd		=  per.odd			,
							U			=  U				,
							n.data		=  n.data			)	
	#	get design matrix
   	design	<- apc.get.design(apc.index,model.design)$design
	##############################
	#	Get point forecasts	for out-of-sample triangle
	linear.predictors.J						<- design.J %*% coefficients.canonical[,1]
	if(distribution.forecasts)
	{
		##############################
		#	transform linear predictor using link to forecasts	
		response.forecast.cell					<- exp(linear.predictors.J)
		trap.response.forecast[index.trap.J]	<- response.forecast.cell
		##############################
		#	Get standard deviations
		#	Use Martinez-Miranda, Nielsen, Nielsen (2005) JRSS-A, Appendix A.2, A.3
		xi.dim	<- nrow(coefficients.canonical)
		tau		<- sum(fitted.values)
		#	MMNN (4.4) pi - multinomial probabilities
		v.pi	<- fitted.values / tau
		#	MMNN (A.2): H.2 insample design
		design.2.average	<- colSums(diag(v.pi) %*% design[,2:xi.dim])
		H.2		<- design[,2:xi.dim] - matrix(data=design.2.average,nrow=n.data,ncol=(xi.dim-1),byrow=T)
		#	MMNN (A.3): information about xi_2 per obs 
	#	i.1.xi.2	<- t(H.2) %*% diag(v.pi) %*% H.2
		#	MMNN (A.12): H.2.J outofsample design
		v.pi.J		<- as.vector(exp(linear.predictors.J) / tau)
		H.2.J		<- design.J[,2:xi.dim] - matrix(data=design.2.average,nrow=n.data.J,ncol=(xi.dim-1),byrow=T)
	#	s2.est.J	<- H.2.J %*% solve(i.1.xi.2,t(H.2.J))
	#	s2.est.J	<- tau * diag(v.pi.J) %*% s2.est.J %*% diag(v.pi.J)
		s2.est.J	<- H.2.J %*% covariance.canonical[2:xi.dim,2:xi.dim] %*% t(H.2.J)
		s2.est.J	<- (tau^2) * diag(v.pi.J) %*% s2.est.J %*% diag(v.pi.J)
		s2.proc.J	<- tau * diag(v.pi.J)
		s2.J		<- s2.proc.J + s2.est.J
		##############################
		#	Function to generate particular forecasts
		function.sum.forecasts <- function(m.X,row.label=NULL)
		{	forecast	<- m.X %*% response.forecast.cell[,1]
			forecast	<- cbind(forecast,sqrt(diag(m.X %*% s2.J	  %*% t(m.X) )))
			forecast	<- cbind(forecast,sqrt(diag(m.X %*% s2.proc.J %*% t(m.X) )))
			forecast	<- cbind(forecast,sqrt(diag(m.X %*% s2.est.J  %*% t(m.X) )))
			colnames(forecast)  <- c("forecast","se","se.proc","se.est")
			if(is.null(row.label)==FALSE)
				rownames(forecast)	<- row.label
			return(forecast)	
		}
		##############################
		#	Generate partial sum forecasts
		#	matrices
		m.partial.sum.age	<- matrix(data=0,nrow=per.forecast.J,ncol=n.data.J)
		m.partial.sum.per	<- matrix(data=0,nrow=per.forecast.J,ncol=n.data.J)
		m.partial.sum.coh	<- matrix(data=0,nrow=per.forecast.J,ncol=n.data.J)
		m.partial.sum.all	<- matrix(data=1,nrow=1,      		 ncol=n.data.J)
		s	<- 0
		for(row in 1:per.forecast.J)
		{	for(age in 1:row)
				m.partial.sum.age[per.forecast.J-row+age,s+age]	<- 1
			for(per in 1:row)
				m.partial.sum.per[per,s+per]		<- 1
			m.partial.sum.coh[row,(s+1):(s+row)]	<- 1	
			s	<- s + row
		}
		#	labels
		label.age	<- c(paste("age_",apc.internal.function.date.2.character(age1+seq(age.max-per.forecast.J,age.max-1)*unit,n.decimal),sep=""))
		label.per	<- c(paste("per_",apc.internal.function.date.2.character(per1+seq(per.max,per.max+per.forecast.J-1)*unit,n.decimal),sep=""))
		label.coh	<- c(paste("coh_",apc.internal.function.date.2.character(coh1+seq(coh.max-per.forecast.J,coh.max-1)*unit,n.decimal),sep=""))
		label.all	<-  "all"                                                                               	
		#	forecasts
		response.forecast.age	<- function.sum.forecasts(m.partial.sum.age,label.age)
		response.forecast.per	<- function.sum.forecasts(m.partial.sum.per,label.per)
		response.forecast.coh	<- function.sum.forecasts(m.partial.sum.coh,label.coh)
		response.forecast.all	<- function.sum.forecasts(m.partial.sum.all,label.all)
		response.forecast.cell	<- function.sum.forecasts(diag(1,n.data.J))
		##############################
		#	Intercept correction
		response.sums		<- apc.data.sums(list(response=response),		data.type="r",apc.index=apc.index)
		m.fitted.values				<- response
		m.fitted.values[index.data]	<- fitted.values
		fitted.values.sums	<- apc.data.sums(list(response=m.fitted.values),data.type="r",apc.index=apc.index)
		intercept.correction.per		<- response.sums$sums.per[per.max] / fitted.values.sums$sums.per[per.max]
		response.forecast.per.ic		<- response.forecast.per
		response.forecast.per.ic[,1]	<- response.forecast.per.ic[,1] * intercept.correction.per
		##############################
		#	Generate partial sum forecasts for sub-groups
		#	by age
		sum.per.by.age.ic <- NULL
		sum.per.by.coh.ic <- NULL
		intercept.correction.per.by.age	<- NULL
		intercept.correction.per.by.coh	<- NULL
		if(is.null(sum.per.by.age)==FALSE)
		{	#	for intercept correction
			data.data.trunc	<- apc.data.list.subset(list(response=response,data.format=data.format,dose=NULL),sum.per.by.age[1]-1,age.max-sum.per.by.age[2],0,0,0,0,apc.index,TRUE)
			data.fit.trunc	<- apc.data.list.subset(list(response=m.fitted.values,data.format=data.format,dose=NULL),sum.per.by.age[1]-1,age.max-sum.per.by.age[2],0,0,0,0,apc.index,TRUE)
			if(sum.per.by.age[1]==sum.per.by.age[2])
			{
				response.sums		<- data.data.trunc$response[length(data.data.trunc$response)]
				fitted.values.sums	<- data.fit.trunc$response[length(data.fit.trunc$response)]
				intercept.correction.per.by.age	<- response.sums/fitted.values.sums
			}
			else
			{	
				response.sums		<- apc.data.sums(data.data.trunc)
				fitted.values.sums	<- apc.data.sums(data.fit.trunc)
				intercept.correction.per.by.age	<- response.sums$sums.per[length(response.sums$sums.per)] / fitted.values.sums$sums.per[length(fitted.values.sums$sums.per)]
			}	
			#	computeting forecast sums	
			sum.length	<- sum.per.by.age[2]-age.max+per.forecast.J
			m.partial.sum.per.by.age	<- matrix(data=0,nrow=sum.length,ncol=n.data.J)
			sum.discard	<- age.max-sum.per.by.age[2]
			s	<- (sum.discard+1)*sum.discard/2
			for(row in (sum.discard+1):per.forecast.J)
			{	per.discard <- max(row-sum.discard-sum.per.by.age[2]+sum.per.by.age[1]-1,0)
				for(per in (1+per.discard):(row-sum.discard))
				 	m.partial.sum.per.by.age[per,s+per]		<- 1
				s	<- s + row
			}
			label.sum	<- c(paste("per_",apc.internal.function.date.2.character(per1+seq(per.max,per.max+sum.length-1)*unit,n.decimal),sep=""))
			sum.per.by.age	<- function.sum.forecasts(m.partial.sum.per.by.age,label.sum)
			sum.per.by.age.ic		<- sum.per.by.age
			sum.per.by.age.ic[,1]	<- sum.per.by.age.ic[,1] * intercept.correction.per.by.age
		}			
		#	by cohort
		if(is.null(sum.per.by.coh)==FALSE)
		{	#	for intercept correction
			data.data.trunc	<- apc.data.list.subset(list(response=response,data.format=data.format,dose=NULL)		,0,0,0,0,sum.per.by.coh[1]-1,coh.max-sum.per.by.coh[2],apc.index,TRUE)
			data.fit.trunc	<- apc.data.list.subset(list(response=m.fitted.values,data.format=data.format,dose=NULL),0,0,0,0,sum.per.by.coh[1]-1,coh.max-sum.per.by.coh[2],apc.index,TRUE)
			if(sum.per.by.coh[1]==sum.per.by.coh[2])
			{
				response.sums		<- data.data.trunc$response[length(data.data.trunc$response)]
				fitted.values.sums	<- data.fit.trunc$response[length(data.fit.trunc$response)]
				intercept.correction.per.by.coh	<- response.sums/fitted.values.sums
			}
			else
			{	
				response.sums		<- apc.data.sums(data.data.trunc)
				fitted.values.sums	<- apc.data.sums(data.fit.trunc)
				intercept.correction.per.by.coh	<- response.sums$sums.per[length(response.sums$sums.per)] / fitted.values.sums$sums.per[length(fitted.values.sums$sums.per)]
			}	
			#	computeting forecast sums	
			sum.length	<- sum.per.by.coh[2]-coh.max+per.forecast.J
			m.partial.sum.per.by.coh	<- matrix(data=0,nrow=sum.length,ncol=n.data.J)
			sum.discard	<- sum.per.by.coh[1]-coh.max+per.forecast.J-1
			s	<- (sum.discard+1)*sum.discard/2
			for(row in 1:(sum.per.by.coh[2]-sum.per.by.coh[1]+1))
			{	for(per in 1:(sum.discard+row))
				 	m.partial.sum.per.by.coh[per,s+per]		<- 1
				s	<- s + sum.discard + row
			}
			label.sum	<- c(paste("per_",apc.internal.function.date.2.character(per1+seq(per.max,per.max+sum.length-1)*unit,n.decimal),sep=""))
			sum.per.by.coh	<- function.sum.forecasts(m.partial.sum.per.by.coh,label.sum)
			sum.per.by.coh.ic		<- sum.per.by.coh
			sum.per.by.coh.ic[,1]	<- sum.per.by.coh.ic[,1] * intercept.correction.per.by.coh
		}			
		##############################
		output.list	<- list(	linear.predictors.forecast		= linear.predictors.J				,
					    		index.trap.J					= index.trap.J						,
					    		trap.response.forecast			= trap.response.forecast			,
					    		response.forecast.cell			= response.forecast.cell			,
					    		response.forecast.age			= response.forecast.age				,
					    		response.forecast.per			= response.forecast.per				,
					    		response.forecast.per.ic		= response.forecast.per.ic			,
					    		response.forecast.coh			= response.forecast.coh				,
					    		response.forecast.all			= response.forecast.all				,
								response.forecast.per.by.age	= sum.per.by.age					,
					    		response.forecast.per.by.age.ic	= sum.per.by.age.ic					,
								response.forecast.per.by.coh	= sum.per.by.coh					,
								response.forecast.per.by.coh.ic	= sum.per.by.coh.ic					,
					    		intercept.correction.per		= intercept.correction.per			,
					    		intercept.correction.per.by.age	= intercept.correction.per.by.age	,
					    		intercept.correction.per.by.coh	= intercept.correction.per.by.coh	)
	}
	if(distribution.forecasts==FALSE)
		return(list(linear.predictors.forecast		= linear.predictors.J	,
					index.trap.J					= index.trap.J			))
	if(distribution.forecasts==TRUE && covariance.output==FALSE)
		return(output.list)
	if(distribution.forecasts==TRUE && covariance.output==TRUE)
		return(c(output.list,
						list(	covariance.proc				= s2.proc.J						,
								covariance.est				= as.vector(exp(linear.predictors.J))
								)))
}	#	apc.forecast.ac

#######################################################################################

apc.forecast.ap <- function(apc.fit,extrapolation.type="I0",suppress.warning=TRUE)
#	BN 31 Jan 2016
#	Computes point forecasts for a model with AP structure.
#	This is done for the triangle which shares age and cohort indices with the data.
#	in
#			apc.fit						List. Output for apc.fit.model
#										Note: apc.fit.model should be run first for a Poisson response model with AC structure so that
#										apc.fit$model.design=="AC"
#										apc.fit$model.family=="poisson.response"
#			extrapolation.type			Character.
#										For the extrapolation future Delta.beta are needed.
#										"I0": Forecasts future beta by average of estimated beta (note, that future Delta.beta invariant to level)
#										"I1": Forecasts future beta by most recent estimated beta
#	out
#			linear.predictors.forecast	Vector. Linear predictors for forecast area
#			index.trap.J				Matrix. age-coh coordinates for vector. Similar structure to index.trap in apc.index
#			trap.response.forecast		Matrix. Includes data and point forecasts. Forecasts in lower right triangle.
{	#	apc.forecast.ap
	########################
	#	get values
	age.max			        <- apc.fit$age.max			   	
	per.max			        <- apc.fit$per.max	
	coh.max			        <- apc.fit$coh.max	
	per.odd			        <- apc.fit$per.odd
	per.zero		        <- apc.fit$per.zero
	U				        <- apc.fit$U
	index.trap		        <- apc.fit$index.trap
	model.design			<- apc.fit$model.design
	linear.predictors		<- apc.fit$linear.predictors
	#	derived values
	per.forecast.J			<- age.max+coh.max-1-per.zero-per.max
	n.data.J				<- (per.forecast.J*(per.forecast.J+1)) %/% 2	# note %/% has higher precedence than *
	##############################
	#	Check Input
	if(isTRUE(model.design=="AP")==FALSE)
		return(cat("ERROR apc.forecast.ap: model.design must be 'AP'\n"))
	if(per.forecast.J==0)
		return(cat("ERROR apc.forecast.ap: forecast area empty since per.forecast=age.max+coh.max-1-per.zero-per.max=0\n"))
	##############################
	#	index for forecast area
	index.trap.J		<- matrix(nrow=n.data.J,ncol=2)
	s	<- 0
	for(row in 1:per.forecast.J)
	{	index.trap.J[(s+1):(s+row),1]	<- seq(age.max-row+1,age.max)	# age
		index.trap.J[(s+1):(s+row),2]	<- coh.max-per.forecast.J+row	# cohort
		s	<- s + row
	}
	colnames(index.trap.J)	<- c("age","coh")
	##############################
	#	Generate forecast trapezoid
	trap.linear.predictors	<- apc.data.list.subset(apc.fit,suppress.warning=suppress.warning)$response
	trap.linear.predictors[index.trap]	<- linear.predictors
	##############################
	#	Extrapolate linear predictor for last diagonal to forecast area
	for(i in 1:per.forecast.J)
		trap.linear.predictors[age.max-per.forecast.J+i,coh.max-seq(0,i-1)] <- trap.linear.predictors[age.max-per.forecast.J+i,coh.max-i]
	##############################
	#	Get design for forecast area
	#	apc index list for forecast area
	apc.index.J	<- list(	age.max		=  age.max			,		
							per.max		=  age.max+coh.max-1,
							coh.max		=  coh.max			,
							index.trap	=  index.trap.J		,
							per.zero	=  per.zero			,
							per.odd		=  per.odd			,
							U			=  U				,
							n.data		=  n.data.J			)
	##############################
	#	Get difference parameters
	id.apc.fit			<- apc.identify(apc.fit)
	coefficients.dif	<- id.apc.fit$coefficients.dif
	index.per.dif		<- id.apc.fit$index.per.dif	
	##############################
	#	Extrapolate period parameters
	#	extrapolation.type == "I0"
	#		form beta_per = beta_1 + sum_{i=2}^per Delta.beta_i , where beta_1 is arbitrary
	#		fit model beta_per = nu
	#		so estimate nu by sum_{per=1}^per.max beta_per
	#		forecast beta_{per.max+h} = nu
	#		thus Delta.beta_{per.max+1} = nu - beta_per.max , which is invariant to beta_1
	#		thus Delpa.beta_{per.max+h} = 0 for h>1
	D.xi.per.extrapolated		<- matrix(data=0,nrow=per.forecast.J,ncol=1)
	if(extrapolation.type=="I0")
	{	D.xi.per	<- coefficients.dif[index.per.dif]
		theta.per	<- cumsum(D.xi.per)
		nu			<- sum(theta.per)/(1+length(D.xi.per))
		D.xi.per.extrapolated[1]	<- nu - theta.per[length(theta.per)]
	}
	#	extrapolation.type == "I1"
	#		form beta_per = beta_1 + sum_{i=2}^per Delta.beta_i , where beta_1 is arbitrary
	#		forecast beta_{per.max+h} = beta_{per.max}
	#		thus Delpa.beta_{per.max+h} = 0 for h>0
	#		thus no change needed
	##############################
	#	partial sum of forecast of Delta.beta
	xi.per.extrapolated	<- cumsum(D.xi.per.extrapolated)
	##############################
	#	update linear predictor
	for(i in 1:per.forecast.J)
		for(j in 1:i)
			trap.linear.predictors[age.max-i+j,coh.max-j+1] <- trap.linear.predictors[age.max-i+j,coh.max-j+1] + xi.per.extrapolated[per.forecast.J-i+1]
	##############################
	return(list(	trap.linear.predictors.forecast	= trap.linear.predictors	,
					index.trap.J				= index.trap.J				,
					D.xi.per.extrapolated		= D.xi.per.extrapolated		))	
}	#	apc.forecast.ap

#######################################################################################


apc.polygon <- function(m.forecast,x.origin=1,
	plot.se=TRUE,plot.se.proc=FALSE,plot.se.est=FALSE,
	unit=1,
	col.line=1,lty.line=1,lwd.line=1,
	q.se=c(2,2,2),
	angle.se=c(45,45,45),
	border.se=c(NA,NA,NA),
	col.se=gray(c(0.50,0.80,0.90)),
	density.se=c(NULL,NULL,NULL),
	lty.se=c(1,1,1))
#	BN 8 Jan 2016
#	Draws point forecasts as line and shaded area to indicate uncertainty
#	This is added to a plot, working like lines() or polygon()
#	In
#        m.forecast	Matrix.  Up to 4 columns.
#        				Column 1: point forecasts.
#								Column 2: forecast standard errors.
#								Column 3: process standard errors.
#								Column 4: estimation standard errors.
#        x.origin		Numerical. x-coordinate for first point forecast. Default: 1.
#        plot.se		Logical. Should forecast   standard errors be plotted? Default: TRUE.
#        plot.se.proc	Logical. Should process    standard errors be plotted? Default: FALSE.
#        plot.se.est	Logical. Should estimation standard errors be plotted? Default: FALSE.
#        unit			Numerical. step length for point forecasts. Default=1.
#        col.line		Point forecasts: Colour of line. Same as col for lines. Default: 1.
#        lty.line		Point forecasts: Type   of line. Same as lty for lines. Default: 1.
#        lwd.line		Point forecasts: Width  of line. Same as lwd for lines. Default: 1.
#        q.se			Vector of length 3. Multiplication factors for standard errors. Default: \code{c(2,2,2)}.}
#        angle.se		Standard error polygon: 3-vector: Angle  of shading.  Same as angle 	for polygon. Default: =c(45,45,45).
#        border.se		Standard error polygon: 3-vector: Border of polygon.  Same as border 	for polygon. Default: =c(NA,NA,NA).
#        col.se			Standard error polygon: 3-vector: Colour of polygon.  Same as col 		for polygon. Default: gray(c(0.50,0.80,0.90)).
#        density.se		Standard error polygon: 3-vector: Density of shading. Same as density 	for polygon. Default: =c(NULL,NULL,NULL).
#        lty.se			Standard error polygon: 3-vector: Type of shading.    Same as lty 		for polygon. Default: =c(1,1,1).
{   #	apc.polygon 
  	x   <- seq(x.origin+unit,by=unit,length=nrow(m.forecast))
  	polygon.forecast <- function(se.choice)
  	{	polygon(c(x,rev(x)),
				c(m.forecast[,1]+q.se[se.choice]*m.forecast[,se.choice+1],
			  rev(m.forecast[,1]-q.se[se.choice]*m.forecast[,se.choice+1])),
				col=col.se[se.choice],
				angle=angle.se[se.choice],
				border=border.se[se.choice],
				density=density.se[se.choice],
				lty=lty.se[se.choice])
  	}
  	if(plot.se)			polygon.forecast(1) 
  	if(plot.se.proc)  	polygon.forecast(2)
  	if(plot.se.est)		polygon.forecast(3)
  	lines(x,m.forecast[,1],col=col.line,lty=lty.line,lwd=lwd.line)        
}	#	apc.polygon
