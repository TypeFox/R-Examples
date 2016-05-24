#######################################################
#	apc package
#	Bent Nielen, 1 April 2015, version 1.0.3
#	functions to fit model
#######################################################
#	Copyright 2014, 2015 Bent Nielsen
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
#	apc.get.design.collinear
#########################################################
apc.get.design.collinear	<- function(apc.index)
#	BN 1 April 2015
#	Constructs a collinear design matrix for an apc model.
#	It includes columns for intercept,
#	for age/period/cohort slopes,  
#	for age/period/cohort double differences.
#	Thus, there are three slopes instead of two.
#	Before use, one has to select which parameters are needed.
#	This should include at either
#	one/two of age/cohort slopes or period slope or no slope.
#	In:		apc.index			List of indices etc.
#	Out:	design.collinear	Matrix.
{	#	apc.get.design.collinear
	########################
	#	get values
	index.trap	<- apc.index$index.trap
	age.max		<- apc.index$age.max
	per.max		<- apc.index$per.max
	coh.max		<- apc.index$coh.max
	per.zero	<- apc.index$per.zero
	per.odd		<- apc.index$per.odd
	U			<- apc.index$U
	n.data		<- apc.index$n.data
	########################
	#	declaring design matrix
	p.design.collinear	<- age.max+per.max+coh.max-2
	  design.collinear	<- matrix(data=0,nrow=n.data,ncol=p.design.collinear)
	########################
	# 	construction of design matrix
	for(row in 1:n.data)
	{
		age	<- index.trap[row,1]		#	age
		coh	<- index.trap[row,2]		#	cohort
		per	<- age+coh-1
		design.collinear[row,1]	<- 1
		design.collinear[row,2]	<- age-U
		design.collinear[row,4]	<- coh-U
		design.collinear[row,3]	<- design.collinear[row,2]+design.collinear[row,4]
		if(age<U)
			design.collinear[row,(4+age):(4+U-1)]										<- seq(1,U-age)
		if(age>U+1)
			design.collinear[row,(4+U):(4+U+age-U-2)]									<- seq(age-U-1,1)
		if(per.odd && per==2*(U-1))
			design.collinear[row,(2+age.max+1)]											<- 1
		if(per>2*U)
			design.collinear[row,(2+age.max+per.odd+1):(2+age.max+per.odd+per-2*U)]		<- seq(per-2*U,1)
		if(coh<U)
			design.collinear[row,(age.max+per.max+coh):(age.max+per.max+U-1)]			<- seq(1,U-coh)
		if(coh>U+1)
			design.collinear[row,(age.max+per.max+U):(age.max+per.max+U+coh-U-2)]		<- seq(coh-U-1,1)
 	}
	return(design.collinear)
}	# 	apc.get.design.collinear

#########################################################
#	apc.get.design
#########################################################
apc.get.design	<- function(apc.index,model.design=NULL)
#	BN 27 Aug 2014
#	Constructs design matrix for an apc model or sub-model thereof.
#	In:		apc.index or apc.fit		List of indices etc.
#			model.design				Character. Indicates which sub-model should be fitted.
#										Default is NULL. Possible choices:
#										"APC","AP","AC","PC","Ad","Pd","Cd","A","P","C","t","tA","tP","tC","1"
#										Argument is not needed when first argument is an apc.fit.
#										If not NULL it will override any information in first argument.
#	Out:	design						Matrix.
#			slopes						Vector. Length 3 of logicals,
#										indicate presence	of age/period/cohort linear slopes
#										at most two slopes can be present
#										if neither age/cohort present then period may be presents,
#										which is the case for model.design "P","tP"
#			difdif						Vector. Length 3 of logicals,
#										indicate presence	of age/period/cohort double differences
{	#	apc.get.design
	##############################
	#	check input
	if(is.null(model.design)==TRUE & is.null(apc.index$model.design)==TRUE)
		return(cat("ERROR apc.get.design: cannot find model.design\n"))
	if(is.null(model.design)==TRUE & is.null(apc.index$model.design)==FALSE)
		model.design	<- apc.index$model.design
	model.design.list		<- c("APC","AP","AC","PC","Ad","Pd","Cd","A","P","C","t","tA","tP","tC","1")
	if(isTRUE(model.design %in% model.design.list)==FALSE)
		return(cat("ERROR apc.get.design: model.design has wrong argument \n"))	
	##############################
	#	get values, that are used
	age.max		<- apc.index$age.max				
	per.max		<- apc.index$per.max				
	coh.max		<- apc.index$coh.max
	##############################
	#	construct indicators for slope parameters and double difference parameters
	#	depending on model.design choice
	#		slopes:		3 vector of logicals for presence of age/period/cohort slope
	#		difdif:		3 vector of logicals for presence of age/period/cohort double differences
	if(model.design=="APC")	{	slopes <- c(1,0,1); difdif <- c(1,1,1);	}
	if(model.design=="AP" )	{	slopes <- c(1,0,1); difdif <- c(1,1,0);	}
	if(model.design=="AC" )	{	slopes <- c(1,0,1); difdif <- c(1,0,1);	}
	if(model.design=="PC" )	{	slopes <- c(1,0,1); difdif <- c(0,1,1);	}
	if(model.design=="Ad" )	{	slopes <- c(1,0,1); difdif <- c(1,0,0);	}
	if(model.design=="Pd" )	{	slopes <- c(1,0,1); difdif <- c(0,1,0);	}
	if(model.design=="Cd" )	{	slopes <- c(1,0,1); difdif <- c(0,0,1);	}
	if(model.design=="A"  )	{	slopes <- c(1,0,0); difdif <- c(1,0,0);	}
	if(model.design=="P"  )	{	slopes <- c(0,1,0); difdif <- c(0,1,0);	}
	if(model.design=="C"  )	{	slopes <- c(0,0,1); difdif <- c(0,0,1);	}
	if(model.design=="t"  )	{	slopes <- c(1,0,1); difdif <- c(0,0,0);	}
	if(model.design=="tA" )	{	slopes <- c(1,0,0); difdif <- c(0,0,0);	}
	if(model.design=="tP" )	{	slopes <- c(0,1,0); difdif <- c(0,0,0);	}
	if(model.design=="tC" )	{	slopes <- c(0,0,1); difdif <- c(0,0,0);	}
	if(model.design=="1"  )	{	slopes <- c(0,0,0); difdif <- c(0,0,0);	}
	##############################
	#	construct index for selecting columns of design matrix
	index.design	<- c(1)
	for(i in 1:3)	if(slopes[i])	index.design	<- c(index.design,i+1)
	if(difdif[1])	index.design <- c(index.design,4+seq(1:(age.max-2)))	
	if(difdif[2])	index.design <- c(index.design,2+age.max+seq(1:(per.max-2)))	
	if(difdif[3])	index.design <- c(index.design,  age.max+per.max+seq(1:(coh.max-2)))
	##############################
	#	get design matrix
	design	<- apc.get.design.collinear(apc.index)[,index.design]
	##############################
	return(list(	design						= design	,
					slopes						= slopes	,
					difdif						= difdif	
			))	
}	#	apc.get.design

#########################################################
#	apc.fit.model
#########################################################
apc.fit.model	<- function(apc.data.list,model.family,model.design,apc.index=NULL)
#	BN  2 feb 2016	Changed: parameter label: date to character changed to allow nice decimal points
#					using apc.internal.function.date.2.character
#	BN 17 mar 2015
#	Function to estimate apc sub-model
#	uses canonical parametrisation as in Nielsen (2013)
#	In:		apc.data.list
#			apc.index	
#			model.family				Character
#										"poisson.response"
#											uses responses only		
#										"od.poisson.response"
#											uses responses only, over-dispersed		
#										"poisson.dose.response"
#											uses doses and responses		
#										"gaussian.response"
#											uses responses only		
#										"gaussian.rates"
#											uses rates=responses/doses only
#										"log.normal.response"
#											takes log of response and fits gaussian model
#			model.design				Character. Indicates which sub-model should be fitted.
#										Possible choices:
#										"APC","AP","AC","PC","Ad","Pd","Cd","A","P","C","t","tA","tP","tC","1"			
#	Out:	fit							List.	Standard output from glm.fit.
#										Note: deviance redefined in Gaussian case,
#										since glm.fit then returns RSS instead of deviance
#			coefficients.canonical		Matrix.	4 columns of coefficients, sdv, t-stat, p-values
#										with rownames for parameters
#			covariance.canonical		Matrix.	Covariance matrix, possibly using mixed parametrisation
#			slopes						Vector. Length 3 of logicals,
#										indicate presence	of age/period/cohort linear slopes
#										at most two slopes can be present
#										if neither age/cohort present then period may be presents,
#										which is the case for model.design "P","tP"
#			difdif						Vector. Length 3 of logicals,
#										indicate presence	of age/period/cohort double differences
#			index.age					Matrix. 1 column. Indices for age    double differences
#			index.per					Matrix. 1 column. Indices for period double differences
#			index.coh					Matrix. 1 column. Indices for cohort double differences
#			dates						Matrix. 1 column. dates for canonical parameters,
#										based on age1,per1,coh1,unit.  Indirectly given by user,
#										through construction of trapezoid format for responses.
#			deviance					Numeric. up to a constant, minus twice the maximized log-likelihood.
#										Where sensible, the constant is chosen so that a saturated model
#										has deviance zero.
#			RSS							Numeric. Residual sum of squares.
#										Only when model.family=="Gaussian.rates"
#			sigma2						Maximum likelihood estimator for variance (divided by n)
#										Only when model.family=="Gaussian.rates"
#			s2							Least squares estimator for variance (divided by df)
#										Only when model.family=="Gaussian.rates"
{	#	apc.fit.model
	##############################
	#	check input
	model.design.list		<- c("APC","AP","AC","PC","Ad","Pd","Cd","A","P","C","t","tA","tP","tC","1")
	model.family.list		<- c("binomial.dose.response","poisson.response","od.poisson.response","poisson.dose.response","gaussian.rates","gaussian.response","log.normal.response")
	model.family.gaussian	<- c("gaussian.rates","gaussian.response","log.normal.response")
	model.family.mixed		<- c("poisson.response","od.poisson.response")
	if(isTRUE(model.design %in% model.design.list)==FALSE)
		return(cat("apc.fit.model error: model.design has wrong argument \n"))	
	if(isTRUE(model.family %in% model.family.list)==FALSE)
		return(cat("apc.fit.model error: model.family has wrong argument \n"))
	######################
	#	get index
	if(is.null(apc.index)==TRUE)
		apc.index	<- apc.get.index(apc.data.list)
	##############################
	#	create indicator for mixed parametrisation
	mixed.par	<- (1-isTRUE(model.design=="1"))*isTRUE(model.family %in% model.family.mixed)
	mixed.par.1	<-    isTRUE(model.design=="1") *isTRUE(model.family %in% model.family.mixed)
	##############################
	#	get values, that are used
	age.max		<- apc.index$age.max				
	per.max		<- apc.index$per.max				
	coh.max		<- apc.index$coh.max
	age1		<- apc.index$age1    
	per1		<- apc.index$per1
	coh1		<- apc.index$coh1    
	unit		<- apc.index$unit
	per.zero	<- apc.index$per.zero
	index.trap	<- apc.index$index.trap
	n.data		<- apc.index$n.data
	v.response	<- apc.index$response[apc.index$index.data]
	v.dose		<- apc.index$dose[    apc.index$index.data]
	n.decimal	<- apc.index$n.decimal
	##############################
	#	get design
	#		design matrix 
	#		slopes:		3 vector of logicals for presence of age/period/cohort slope
	#		difdif:		3 vector of logicals for presence of age/period/cohort double differences
	get.design	<- apc.get.design(apc.index,model.design)
	design		<- get.design$design
	slopes		<- get.design$slopes
	difdif		<- get.design$difdif
	##############################
	#	REGRESSION
	#	Binomial/logistic regression with logit link
	if(model.family=="binomial.dose.response")
		fit	<- glm.fit(design,cbind(v.response,v.dose-v.response),family=binomial(link="logit"))
	#	Poisson regression for response only and with log link
	if(model.family=="poisson.response")
		fit	<- glm.fit(design,v.response,family=poisson(link="log"))	
	if(model.family=="od.poisson.response")
		fit	<- glm.fit(design,v.response,family=poisson(link="log"))	
	#	Poisson regression for dose-response and with log link
	if(model.family=="poisson.dose.response")
		fit	<- glm.fit(design,v.response,family=poisson(link="log"),offset=log(v.dose))
	#	Gaussian regression for response only and with identity link (Least Squares)
	if(model.family=="gaussian.response")
		fit	<- glm.fit(design,v.response,family=gaussian(link="identity"))
	#	Gaussian regression for response only and with identity link (Least Squares)
	if(model.family=="gaussian.rates")
		fit	<- glm.fit(design,v.response/v.dose,family=gaussian(link="identity"))		
	#	Gaussian regression for log(response) and with identity link (Least Squares)
	if(model.family=="log.normal.response")
		fit	<- glm.fit(design,log(v.response),family=gaussian(link="identity"))
	##############################
	#	construct for indices for double difference parameters  
	index.age	<- NULL
	index.per	<- NULL
	index.coh	<- NULL
	start		<- 1+sum(slopes)
	if(difdif[1])	{	index.age	<- start+seq(1,age.max-2);	start	<- start+age.max-2	}
	if(difdif[2])	{	index.per	<- start+seq(1,per.max-2);	start	<- start+per.max-2	}
	if(difdif[3])	{	index.coh	<- start+seq(1,coh.max-2);	start	<- start+coh.max-2	}
	xi.dim		<- start
	##############################
	#	construct dates for for double difference parameters
	dates		<- matrix(data=NA,nrow=xi.dim,ncol=1)			
	if(difdif[1])	dates[index.age,1]	<- age1+seq(2,age.max-1)*unit	
	if(difdif[2])	dates[index.per,1]	<- per1+seq(2,per.max-1)*unit
	if(difdif[3])	dates[index.coh,1]	<- coh1+seq(2,coh.max-1)*unit
	##############################
	#	construct row names for canonical parameter
	names	<- c("level")
	if(slopes[1])	names	<- c(names,"age slope")
	if(slopes[2])	names	<- c(names,"period slope")
	if(slopes[3])	names	<- c(names,"cohort slope")
	if(difdif[1])
		for(i in 1:(age.max-2))
			names	<- c(names,paste("DD_age_"   ,apc.internal.function.date.2.character((dates[index.age,1])[i],n.decimal),sep=""))
			if(difdif[2])
		for(i in 1:(per.max-2))
			names	<- c(names,paste("DD_period_",apc.internal.function.date.2.character((dates[index.per,1])[i],n.decimal),sep=""))
	if(difdif[3])
		for(i in 1:(coh.max-2))
			names	<- c(names,paste("DD_cohort_",apc.internal.function.date.2.character((dates[index.coh,1])[i],n.decimal),sep=""))
	##############################
	#	Get coefficients and covariance
	coefficients.canonical	<- summary.glm(fit)$coefficients
	rownames(coefficients.canonical)	<- names
	covariance.canonical	<- summary.glm(fit)$cov.scaled
	#	Need to condition in mixed parametrisation
	if(mixed.par)
	{
		c22	<- covariance.canonical[2:xi.dim,2:xi.dim]
		c21	<- covariance.canonical[2:xi.dim,1]
		c11	<- covariance.canonical[1,1]
		covariance.canonical[2:xi.dim,2:xi.dim]	<- c22 - c21 %o% c21 / c11
	}
	if(mixed.par | mixed.par.1)
	{
		covariance.canonical[1,]	<- 0
		covariance.canonical[,1]	<- 0		
	}	
	##############################
	#	get standard errors 
	coefficients.canonical[,2]	<- sqrt(diag(covariance.canonical))
	if(mixed.par | mixed.par.1)		# mixed parametrisation
		coefficients.canonical[1,2]	<- NA
	##############################
	#	get t-statistics
	coefficients.canonical[,3]	<- coefficients.canonical[,1] 	/ coefficients.canonical[,2]
	##############################
	#	get p-values
	coefficients.canonical[,4]	<- 2*pnorm(abs(coefficients.canonical[  ,3]),lower.tail=FALSE)
	##############################
	#	get deviance, RSS and variance estimate in Gaussian case
	#	note that glm.fit returns RSS instead of deviance in Gaussian case
	RSS		<- NULL
	sigma2	<- NULL
	s2		<- NULL
	if(isTRUE(model.family %in% model.family.gaussian))
	{
		RSS				<- fit$deviance
		sigma2			<- RSS/n.data
		s2				<- RSS/fit$df.residual
		fit$deviance	<- n.data*(1+log(2*pi)+log(sigma2))
	}	
	##############################
	return(c(fit,apc.index,
				list(	model.family				= model.family					,
						model.design				= model.design					,
						coefficients.canonical		= coefficients.canonical		,
						covariance.canonical		= covariance.canonical			,
						slopes						= slopes						,
						difdif						= difdif						,
						index.age					= index.age 					,
						index.per					= index.per 					,
						index.coh					= index.coh 					,
						dates						= dates							,
						RSS							= RSS							,
						sigma2						= sigma2						,
						s2							= s2							
			)))	
}	#	apc.fit.model

#########################################################
#	apc.fit.table
#########################################################
apc.fit.table	<- function(apc.data.list,model.family,apc.index=NULL)
#	BN 17 mar 2015
{	#	apc.fit.table
	######################
	#	model families
	model.family.list		<- c("binomial.dose.response","poisson.response","od.poisson.response","poisson.dose.response","gaussian.rates","gaussian.response","log.normal.response")
	model.family.gaussian	<- c("gaussian.rates","gaussian.response","log.normal.response")
	model.family.od			<- c("od.poisson.response")
	######################
	#	check input
	if(isTRUE(model.family %in% model.family.list)==FALSE)
		return(cat("apc.fit.table error: model.family has wrong argument \n"))
	######################
	#	get index
	if(is.null(apc.index)==TRUE)
		apc.index	<- apc.get.index(apc.data.list)
	######################
	#	Function to get one line of table from two fits
	fit.tab.line.glm	<- function(fit.U,fit.R,gaussian)
	#	BN 20 Sep 2013
	{
		dev.U	<- fit.U$deviance
		dev.R	<- fit.R$deviance
		df.U	<- fit.U$df.residual
		df.R	<- fit.R$df.residual
		LR	<- dev.R-dev.U
		df	<- df.R-df.U
		aic	<- fit.R$aic
		if(gaussian)		
			return(round(c(dev.R,df.R,LR,df,pchisq(LR,df,lower.tail=FALSE),aic),digits=3))
		else	
			return(round(c(dev.R,df.R,pchisq(dev.R,df.R,lower.tail=FALSE),LR,df,pchisq(LR,df,lower.tail=FALSE),aic),digits=3))
	}
	######################	
	model.design.list	<- c("APC","AP","AC","PC","Ad","Pd","Cd","A","P","C","t","tA","tP","tC","1")
	######################	
	#	number of columns
															ncol <- 7
		if(isTRUE(model.family %in% model.family.gaussian))	ncol <- 6
		if(isTRUE(model.family %in% model.family.od))		ncol <- 9
	######################
	#	declare table
	fit.tab		<- matrix(nrow=length(model.design.list),ncol=ncol,data=NA)
	#	unrestricted apc model
	fit.apc	<- apc.fit.model(apc.data.list,model.family,model.design="APC",apc.index)	
	#	model list
	for(i in 1:length(model.design.list))
	{
		fit.sub			<- apc.fit.model(apc.data.list,model.family,model.design=model.design.list[i],apc.index)
		if(isTRUE(model.family %in% model.family.gaussian))
			fit.tab[i,1:6]		<- fit.tab.line.glm(fit.apc, fit.sub,1)
		else	
			fit.tab[i,1:7]		<- fit.tab.line.glm(fit.apc, fit.sub,0)
	}
	#	OD F-test
	if(isTRUE(model.family %in% model.family.od))
	for(i in 2:length(model.design.list))
	{	fit.tab[i,8]	= (fit.tab[i,4]/fit.tab[i,5])/(fit.tab[1,1]/fit.tab[1,2])
		fit.tab[i,9]	= round(pf(fit.tab[i,8],fit.tab[i,5],fit.tab[1,2],lower.tail=FALSE),digits=3)
	}
	#	insert NAs
		insert			<- c(4,5,6)
	if(isTRUE(model.family %in% model.family.gaussian))
		insert			<- c(3,4,5)
	fit.tab[1,insert] <- NA
	#	row names
	rownames(fit.tab)	<- model.design.list
	#	column names
		colnames		<- c("-2logL","df.residual","prob(>chi_sq)","LR.vs.APC","df.vs.APC","prob(>chi_sq)","aic")
	if(isTRUE(model.family %in% model.family.gaussian))
		colnames		<- c("-2logL","df.residual","LR.vs.APC","df.vs.APC","prob(>chi_sq)","aic")
	if(isTRUE(model.family %in% model.family.od))
		colnames		<- c("-2logL","df.residual","prob(>chi_sq)","LR.vs.APC","df.vs.APC","prob(>chi_sq)","aic","F","prob(>F)")	
	colnames(fit.tab)	<- colnames		
	######################
	return(fit.tab)
}	#	apc.fit.table