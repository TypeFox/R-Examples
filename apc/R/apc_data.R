#######################################################
#	apc package
#	Bent Nielsen, 6 Feb 2016, version 1.2
#	Data list and Data examples
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

apc.data.list	<- function(response, data.format, dose=NULL, age1=1, per1=1, coh1=1, unit=1, per.zero=NULL, per.max=NULL, time.adjust=0, label=NULL, n.decimal=NULL)
#	BN 1 Feb 2016: 	Added argument n.decimal
#	BN 8 Sep 2015:	Added argument label
#	BN 24 apr 2015
#	This function constructs list of apc.data.list type.
#	This gives the user a single focus for entering information about the data.
#	Only response and data.format are obligatory input.
#	in:		response		matrix of responses
#			dose			NULL or matrix of dose.
#			data.format		character indicating format of data.matrix
#									"AP"		has    age/period as increasing row/column index
#									"AC"		has    age/cohort as increasing row/column index
#									"CA"		has cohort/age    as increasing row/column index
#									"CL"		has cohort/age 	  as increasing row/column index, triangular
#									"CP"		has cohort/period as increasing row/column index
#									"PA"		has period/age    as increasing row/column index
#									"PC"		has period/cohort as increasing row/column index
#									"trapezoid"	has    age/period as increasing row/column index,
#													period-diagonals are NA for period <= per.zero and >per.zero+per.max 
#			age1			smallest age    index (not used for data.format="CP", "PC")
#			per1			smallest period index (not used for data.format="AC","CA","CL","CL.vector.by.row","trapezoid")
#			coh1			smallest cohort index (not used for data.format="AP","PA")
#			unit			time units for indices
#  			per.zero		Only used for data.format="trapezoid".
#									Entries in upper triangle with period <= per.zero are interpreted as NA.
#  			per.max 		Only used for data.format="trapezoid".
# 	 								Entries in lower triangle with period > per.zero+per.max are interpreted as NA.
#			time.adjust		Only two of age1, per1, coh1 are used.
#							The third is computed according to the formula
#							age1+coh1=per1+time.adjust
#			label			character.
#							particularly useful with multiple data sets
#	out		list including all 8 arguments.
{	#	apc.data.list
	##############################
	#	check obligatory input
	data.format.list		<- c("AP","AC","CA","CL","CP","PA","PC","trap","trapezoid")
	data.format.list.matrix	<- c("AP","AC","CA","CL","CP","PA","PC","trap","trapezoid")
	if(isTRUE(data.format %in% data.format.list)==FALSE)
		return(cat("apc.error: model.family has wrong argument \n"))
	if(is.matrix(response)==FALSE)
		return(cat("apc.error: response is not a matrix \n"))
	#	check "CL" input
	if(data.format=="CL")
	{
		if(ncol(response) != nrow(response))	return(cat("apc.error: Response matrix is not square \n"))
		k	<- nrow(response)
		for(age in 2:k)
			for(coh in (k+2-age):k)
				if(is.na(response[coh,age])==FALSE) return(cat("apc.error: Lower triangle of response matrix should be NA \n"))
	}		
	##############################
	if(data.format=="trap")			data.format	<- "trapezoid"
	if(data.format=="trapezoid")
	{
		if(is.null(per.zero))	per.zero	<-0;
		if(is.null(per.max))	per.max	<-nrow(response)+ncol(response)-1-per.zero;
	}
	else
	{	per.zero	<- NULL;
		per.max	<- NULL;
	}
	##############################
	#	row and column names
	function.get.dim.names <- function(m,x1,x.unit,y1,y.unit)
	{	#	function.get.dim.names
		function.one.set.of.names <- function(mm,dim1,dim.unit)
		{	#	function.one.set.of.names
			dim.length	<- nrow(mm)
#			if(is.integer(dim.unit) && dim.unit>1)
#				dim.names	<- paste(as.character(seq(from=dim1,length=dim.length,by=dim.unit)),
#							 		 "-",
#							 		 as.character(seq(from=dim1,length=dim.length,by=dim.unit)+dim.unit-1),
#									 sep="")
#			else									 
				dim.names	<- as.character(seq(from=dim1,length=dim.length,by=dim.unit))
			return(dim.names)						 
		}	#	function.one.set.of.names
		if(is.null(rownames(m)))
			rownames(m) <- function.one.set.of.names(m,x1,x.unit)
		if(is.null(colnames(m)))
			colnames(m) <- function.one.set.of.names(t(m),y1,y.unit)
		return(m)
	}	#	function.get.dim.names
	if(data.format %in% c("AC","trap","trapezoid"))
		{	x1 <- age1; y1 <- coh1	}
	if(data.format %in% c("CA","CL"))
		{	x1 <- coh1; y1 <- age1	}	
	if(data.format %in% c("AP"))
		{	x1 <- age1; y1 <- per1	}	
	if(data.format %in% c("PA"))
		{	x1 <- per1; y1 <- age1	}	
	if(data.format %in% c("CP"))
		{	x1 <- coh1; y1 <- per1	}	
	if(data.format %in% c("PC"))
		{	x1 <- per1; y1 <- coh1	}	
	response 	<- function.get.dim.names(response,x1,unit,y1,unit)
	if(is.null(dose)==FALSE)
		dose	<- function.get.dim.names(dose,	x1,unit,y1,unit)
	##############################
	#	n decimal
	if(is.null(n.decimal))
		if(unit<1 && unit>=1/20)	n.decimal	<- 2
	##############################
	return(list(response	=response	,
				dose		=dose		,
				data.format	=data.format,
				age1		=age1		,
				per1		=per1		,
				coh1		=coh1		,
				unit		=unit		,
				per.zero	=per.zero	,
				per.max		=per.max	,
				time.adjust	=time.adjust,
				label		=label		,
				n.decimal	=n.decimal	))
}	#	apc.data.list

vector.2.triangle	<- function(v,k)
#	BN 7 Feb 2015
#	function to organise a vector as a triangle.
#	useful for reserving data
#	in:		v		vector. Length k*(k+1)/2
#			k		integer. Dimension	
#	out:	m		matrix with "CL" format Dimension kxk.
#					Upper left triangle filled by v, row by row.
#					Remaining entries NA
{	#	vector.2.triangle	
	##############################
	#	Check input
	if(is.vector(v)==FALSE)		return(cat("vector.2.triangle: v is not a vector \n"))
	if(length(v) != k*(k+1)/2)	return(cat("vector.2.triangle: Length of v does not match k\n"))
	##############################
	#	turn into matrix
	m	<- matrix(nrow=k,ncol=k,data=NA)
	i	<- 0
	for(coh in 1:k)
	{
		m[coh,1:(k+1-coh)]	<- v[(i+1):(i+k+1-coh)]
		i	<- i+k+1-coh
	}
	return(m)
}	#	vector.2.triangle	

#############################################################################################
#	SPECIFIC DATA SETS
#############################################################################################

###############################
#	JAPANESE BREAST CANCER DATA
###############################
data.Japanese.breast.cancer	<- function()
#	BN, 24 apr 2015	(17 oct 2013)
#	An example with A,P,C effects
#
#	Taken from table I of
#	Clayton, D. and Schifflers, E. (1987b)
#	Models for temperoral variation in cancer rates. II: age-period-cohort models.
#	Statistics in Medicine 6, 469-481.
#
#	Original Table caption:
#	age-specific mortality rates (per 100,000 person-years observation) of breast cancer in Japan,
#	during the period 1955-1979. Numbers of cases on which rates are based are in parentheses
#	(source: WHO mortality data base).  
{	#	data.Japanese.breast.cancer
v.rates		<- c( 0.44, 0.38, 0.46, 0.55, 0.68,
			 	  1.69, 1.69, 1.75, 2.31, 2.52,
				  4.01, 3.90, 4.11, 4.44, 4.80,
				  6.59, 6.57, 6.81, 7.79, 8.27,
				  8.51, 9.61, 9.96,11.68,12.51,
				 10.49,10.80,12.36,14.59,16.56,
				 11.36,11.51,12.98,14.97,17.79,
				 12.03,10.67,12.67,14.46,16.42,
				 12.55,12.03,12.10,13.81,16.46,
				 15.81,13.87,12.65,14.00,15.60,
				 17.97,15.62,15.83,15.71,16.52)
v.cases		<- c(   88,   78,  101,  127,  179,
				   299,  330,  363,  509,  588,
				   596,  680,  798,  923, 1056,
				   874,  962, 1171, 1497, 1716,
				  1022, 1247, 1429, 1987, 2398,
				  1035, 1258, 1560, 2079, 2794,
				   970, 1087, 1446, 1828, 2465,
				   820,  861, 1126, 1549, 1962,
				   678,  738,  878, 1140, 1683,
				   640,  628,  656,  900, 1162,
				   497,  463,  536,  644,  865)				 
col.names	<- paste(as.character(seq(from=1955,length=5,by=5)),"-",
			 		 as.character(seq(from=1955,length=5,by=5)+4),sep="")
row.names	<- paste(as.character(seq(from=25  ,length=11,by=5)),"-",
			 		 as.character(seq(from=25  ,length=11,by=5)+4),sep="")
rates	<- matrix(data=v.rates,nrow=11, ncol=5,byrow=TRUE,dimnames=list(row.names,col.names))
cases	<- matrix(data=v.cases,nrow=11, ncol=5,byrow=TRUE,dimnames=list(row.names,col.names))
return(apc.data.list(
			response	=cases			,
			dose		=cases/rates	,
			data.format	="AP"			,
			age1		=25				,
			per1		=1955			,
			unit		=5				,
			label		="Japanese breast cancer"))
}	#	data.Japanese.breast.cancer

##################################
#	ITALIAN BLADDER CANCER DATA
##################################
data.Italian.bladder.cancer	<- function()
#	BN, 24 apr 2015 (17 oct 2013)
#	An example with A,C effects
#
#	Taken from table IV of
#	Clayton, D. and Schifflers, E. (1987a)
#	Models for temperoral variation in cancer rates. I: age-period and age-cohort models.
#	Statistics in Medicine 6, 449-467.
#
#	Original Table caption:
#	age-specific incidence rates (per 100,000 person-years observation) of bladder cancer in
#	Italian males during the period 1955-1979. Numerators are in parentheses
#	(source of data: WHO mortality database).  
{	#	data.Italian.bladder.cancer
v.rates		<- c( 0.03, 0.03, 0.01, 0.04,  0.12,
				  0.17, 0.18, 0.12, 0.08,  0.09,
				  0.32, 0.31, 0.35, 0.42,  0.32,
				  1.04, 1.05, 0.91, 1.04,  1.27,
				  2.86, 2.52, 2.61, 3.04,  3.16,
				  6.64, 7.03, 6.43, 6.46,  8.47,
				 12.71,13.39,14.59,14.64, 16.38,
				 20.11,23.98,26.69,27.55, 28.53,
				 24.40,33.16,42.12,47.77, 50.37,
				 32.81,42.31,52.87,66.01, 74.64,
				 45.54,47.94,62.05,84.65,104.21)
v.cases		<- c(   3,   3,   1,   4,  12,
  				   16,  17,  11,   8,   8,
				   24,  29,  33,  39,  30,
				   79,  76,  82,  95, 115,
				  234, 185, 183, 267, 285,
				  458, 552, 450, 431, 723,
				  720, 867,1069, 974,1004,
				  890,1230,1550,1840,1811,
				  891,1266,1829,2395,3028,
				  920,1243,1584,2292,3176,
				  831, 937,1285,1787,2659)
col.names	<- paste(as.character(seq(from=1955,length=5,by=5)),"-",
			 		 as.character(seq(from=1955,length=5,by=5)+4),sep="")
row.names	<- paste(as.character(seq(from=25  ,length=11,by=5)),"-",
			 		 as.character(seq(from=25  ,length=11,by=5)+4),sep="")
rates	<- matrix(data=v.rates,nrow=11, ncol=5,byrow=TRUE,dimnames=list(row.names,col.names))
cases	<- matrix(data=v.cases,nrow=11, ncol=5,byrow=TRUE,dimnames=list(row.names,col.names))
return(apc.data.list(
			response	=cases			,
			dose		=cases/rates	,
			data.format	="AP"			,
			age1		=25				,
			per1		=1955			,
			unit		=5				,
			label		="Italian bladder cancer"))
}	#	data.Italian.bladder.cancer

##################################
#	BELGIAN LUNG CANCER DATA
##################################
data.Belgian.lung.cancer	<- function(unbalanced=FALSE)
#	BN, 17 oct 2013
#	An example with A,drift effects
#
#	Taken from table VIII of
#	Clayton, D. and Schifflers, E. (1987a)
#	Models for temperoral variation in cancer rates. I: age-period and age-cohort models.
#	Statistics in Medicine 6, 449-467.
#
#	Original Table caption:
#	age-specific mortality rates (per 100,000 person-years observation) of lung cancer in
#	Belgian females during the period 1955-1978. Numerators are shown in parentheses
#	(source of data: WHO mortality database).
#
#	NOTE	The data are unbalanced since the last column only covers 4 years.  This is not used.
#	In:		unbalanced		logical.  If true unbalanced version includind last column
{	#	data.Belgian.lung.cancer
v.rates		<- c( 0.19, 0.13, 0.50, 0.19, 0.70,
				  0.66, 0.98, 0.72, 0.71, 0.57,
				  0.78, 1.32, 1.47, 1.64, 1.32,
				  2.67, 3.16, 2.53, 3.38, 3.93,
				  4.84, 5.60, 4.93, 6.05, 6.83,
				  6.60, 8.50, 7.65,10.59,10.42,
				 10.36,12.00,12.68,14.34,17.95,
				 14.76,16.37,18.00,17.60,23.91,
				 20.53,22.60,24.90,24.33,32.70,
				 26.24,27.70,30.47,36.94,38.47,
				 33.47,33.61,36.77,43.69,45.20)
v.cases		<- c(  3,  2,  7,  3, 10,
				  11, 16, 11, 10,  7,
				  11, 22, 24, 25, 15,
				  36, 44, 42, 53, 48,
				  77, 74, 68, 99, 88,
				 106,131, 99,142,134,
				 157,184,189,180,177,
				 193,232,262,249,239,
				 219,267,323,325,343,
				 223,250,308,412,358,
				 198,214,253,338,312)
col.names	<- c("1955-1959","1960-1964","1965-1969","1970-1974","1975-1978")
row.names	<- paste(as.character(seq(from=25  ,length=11,by=5)),"-",
			 		 as.character(seq(from=25  ,length=11,by=5)+4),sep="")
rates	<- matrix(data=v.rates,nrow=11, ncol=5,byrow=TRUE,dimnames=list(row.names,col.names))
cases	<- matrix(data=v.cases,nrow=11, ncol=5,byrow=TRUE,dimnames=list(row.names,col.names))
if(unbalanced==FALSE)
	return(apc.data.list(
			response	=cases[,(1:4)]					,
			dose		=cases[,(1:4)]/rates[,(1:4)]	,
			data.format	="AP"							,
			age1		=25								,
			per1		=1955							,
			unit		=5								))
if(unbalanced==TRUE)
	return(apc.data.list(
			response	=cases			,
			dose		=cases/rates	,
			data.format	="AP"			,
			unit		=5				,
			label		="Belgian lung cancer"))
}	#	data.Belgian.lung.cancer

###################################
##	PROSTATE CANCER FOR NONWHITES IN THE US
###################################
data.US.prostate.cancer	<- function()
##	BN, 28 apr 2015
##	An example with over-dispersion
##
##	Taken from table 2 of
##	Holford, T.R. (1983)
##	The estimation of age, period and cohort effects for vital rates.
##	Biometrics 39, 311-324.
##
##	Original Table caption:
##	Number of prostate cancer deathrs and midperiod population for nonwhites in the
##	U.S. by age and period
##	Sources:
##	Cancer deaths: National Center for Health Statistics, 1937-1973
##	Population 1935-60: Grove and Hetzel, 1968
##	Population 1960-69: Bureau of the Census, 1974		
##	Population measured in 1000s
##
{	#	data.US.prostate.cancer
v.deaths	<- c( 177, 271, 312, 382, 321, 305, 308,
				  262, 350, 552, 620, 714, 649, 738,
				  360, 479, 644, 949, 932,1292,1327,
				  409, 544, 812,1150,1668,1958,2153,
				  328, 509, 763,1097,1593,2039,2433,
				  222, 359, 584, 845,1192,1638,2068,
				  108, 178, 285, 475, 742, 992,1374)

v.population<- c( 301, 317, 353, 395, 426, 473, 498,
				  212, 248, 279, 301, 358, 411, 443,
				  159, 194, 222, 222, 258, 304, 341,
				  132, 144, 169, 210, 230, 264, 297,
				   76,  94, 110, 125, 149, 180, 197,
				   37,  47,  59,  71,  91, 108, 118,
				   19,  22,  32,  39,  44,  56,  66)
col.names	<- paste(as.character(seq(from=1935,length=7,by=5)),"-",
			 		 as.character(seq(from=1935,length=7,by=5)+4),sep="")
row.names	<- paste(as.character(seq(from=50  ,length=7,by=5)),"-",
			 		 as.character(seq(from=50  ,length=7,by=5)+4),sep="")
response	<- matrix(data=v.deaths		,nrow=7, ncol=7,byrow=TRUE,dimnames=list(row.names,col.names))
dose		<- matrix(data=v.population	,nrow=7, ncol=7,byrow=TRUE,dimnames=list(row.names,col.names))
return(apc.data.list(
			response	=response		,
			dose		=dose			,
			data.format	="AP"			,
			age1		=50				,
			per1		=1935			,
			unit		=5				,
			label		="US prostate cancer"))
}	#	data.US.prostate.cancer

##################################
#	UK Asbestos data
##################################
data.asbestos	<- function(all.age.groups=FALSE)
#	BN, 17 oct 2013
#
#	Taken from
#	Martinez Miranda, Nielsen and Nielsen (2013)
#	Inference and forecasting in the age-period-cohort model with unknown exposure with
#	an application to mesothelioma mortality.
#	To appear in Journal of the Royal Statistical Society series A
#
#	
#	update of data from the Health and Safety Executive
{	#	data.asbestos
v.cases	<-c(0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,1	,1	,0	,2	,0	,1	,2	,0	,0	,1	,0	,0	,2	,1	,1	,1	,1	,4	,1	,1	,4	,5	,3	,5	,3	,3	,6	,3	,2	,3	,4	,1	,4	,1	,0	,2	,1	,0	,1	,0	,0	,0	,0	,2	,0	,0	,0	,1	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,1	,1	,0	,0	,1	,0	,0	,0	,0	,0	,1	,2	,0	,3	,0	,1	,3	,4	,1	,1	,2	,6	,1	,1	,3	,3	,5	,3	,4	,1	,5	,3	,8	,3	,4	,4	,5	,3	,1	,3	,2	,2	,3	,1	,0	,3	,1	,4	,2	,0	,1	,1	,3	,1	,0	,0	,1	,0	,1	,1	,0	,0	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,1	,3	,2	,2	,3	,0	,1	,1	,1	,0	,3	,1	,4	,6	,8	,3	,6	,3	,10	,6	,7	,6	,4	,2	,6	,5	,8	,4	,1	,0	,2	,1	,3	,1	,1	,1	,2	,1	,0	,0	,0	,0	,1	,0	,1	,0	,0	,0	,1	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,2	,1	,0	,1	,0	,0	,1	,1	,2	,1	,0	,1	,1	,4	,2	,1	,2	,5	,4	,6	,5	,10	,3	,4	,11	,10	,5	,9	,1	,5	,4	,7	,6	,3	,2	,4	,6	,5	,0	,1	,0	,0	,1	,2	,2	,1	,0	,1	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,2	,0	,0	,1	,1	,1	,1	,1	,3	,0	,1	,0	,3	,3	,4	,3	,2	,4	,3	,5	,5	,1	,1	,10	,7	,4	,7	,5	,2	,5	,13	,1	,5	,3	,4	,0	,6	,5	,1	,4	,2	,2	,1	,1	,1	,2	,1	,0	,0	,1	,1	,1	,0	,0	,0	,0	,0	,0	,0,
			0	,0	,1	,0	,0	,0	,0	,0	,1	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,2	,1	,2	,1	,6	,1	,2	,5	,2	,4	,4	,2	,3	,4	,6	,9	,7	,8	,8	,3	,6	,6	,5	,2	,6	,7	,4	,10	,5	,5	,3	,5	,6	,2	,1	,1	,2	,0	,3	,0	,1	,1	,0	,1	,1	,0	,2	,0	,0	,0	,0	,0	,0	,1,
			1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,2	,0	,1	,0	,1	,0	,0	,1	,1	,1	,1	,1	,1	,0	,4	,4	,1	,5	,5	,6	,1	,6	,5	,2	,6	,1	,5	,8	,5	,9	,9	,6	,7	,8	,5	,3	,7	,9	,7	,4	,8	,2	,5	,4	,2	,1	,4	,2	,0	,1	,0	,1	,1	,1	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,2	,0	,0	,0	,2	,0	,1	,3	,4	,3	,5	,2	,6	,2	,5	,3	,4	,4	,11	,3	,5	,10	,10	,3	,6	,11	,7	,8	,6	,6	,4	,9	,10	,7	,5	,2	,3	,2	,0	,4	,0	,0	,2	,2	,0	,2	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,1	,1	,1	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,2	,1	,2	,2	,1	,2	,6	,2	,3	,12	,7	,5	,3	,4	,3	,4	,3	,8	,8	,6	,11	,11	,9	,11	,11	,4	,6	,10	,5	,7	,6	,9	,3	,3	,3	,3	,5	,0	,4	,2	,3	,1	,1	,0	,0	,0	,1	,0	,1	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,1	,0	,1	,1	,2	,0	,0	,2	,2	,2	,0	,1	,0	,3	,3	,3	,7	,5	,4	,5	,9	,5	,8	,9	,5	,7	,5	,14	,13	,5	,11	,9	,7	,10	,8	,9	,9	,12	,8	,2	,11	,7	,7	,3	,0	,4	,3	,3	,1	,2	,3	,1	,0	,2	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,1	,0	,0	,2	,0	,3	,2	,2	,1	,3	,3	,5	,6	,1	,5	,7	,5	,6	,5	,6	,5	,11	,9	,4	,10	,4	,9	,9	,9	,14	,13	,10	,7	,6	,8	,10	,10	,8	,7	,7	,9	,8	,2	,4	,2	,2	,1	,3	,2	,1	,1	,1	,1	,0	,0	,0	,0	,1	,0	,0	,0,
			0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,1	,1	,0	,2	,1	,0	,0	,1	,1	,2	,2	,2	,2	,1	,5	,3	,7	,5	,5	,9	,8	,9	,13	,11	,9	,8	,8	,12	,11	,9	,12	,6	,23	,5	,17	,11	,8	,4	,5	,8	,13	,12	,12	,9	,8	,3	,5	,4	,6	,3	,1	,0	,1	,1	,1	,0	,0	,2	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,1	,0	,1	,1	,1	,2	,1	,0	,1	,3	,3	,1	,0	,5	,0	,4	,8	,4	,7	,10	,10	,9	,9	,12	,11	,10	,10	,8	,8	,6	,8	,14	,10	,13	,13	,15	,15	,10	,13	,15	,8	,12	,8	,11	,6	,6	,6	,3	,1	,2	,2	,2	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,2	,1	,0	,2	,0	,0	,3	,1	,2	,1	,1	,6	,2	,7	,3	,5	,5	,6	,11	,11	,13	,8	,8	,13	,12	,17	,9	,15	,8	,6	,10	,13	,17	,16	,14	,12	,11	,10	,9	,12	,8	,4	,9	,5	,7	,7	,4	,0	,1	,1	,2	,1	,3	,1	,0	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,1	,1	,1	,1	,2	,4	,3	,3	,1	,1	,4	,4	,3	,6	,3	,4	,10	,3	,9	,10	,17	,12	,13	,14	,18	,17	,11	,14	,18	,12	,12	,16	,14	,12	,12	,11	,12	,5	,14	,9	,7	,11	,12	,3	,7	,7	,5	,3	,1	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,1	,2	,0	,1	,2	,1	,4	,3	,3	,4	,4	,2	,3	,5	,5	,5	,3	,4	,4	,11	,10	,7	,14	,5	,18	,13	,15	,12	,22	,11	,13	,10	,15	,21	,12	,14	,14	,16	,22	,15	,6	,14	,6	,11	,8	,5	,4	,2	,1	,3	,3	,2	,1	,1	,0	,1	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,3	,0	,1	,2	,1	,5	,6	,3	,6	,4	,4	,5	,10	,7	,9	,13	,10	,12	,16	,14	,21	,21	,18	,12	,16	,11	,11	,5	,20	,24	,14	,21	,11	,15	,20	,14	,17	,11	,9	,7	,7	,9	,6	,12	,2	,3	,1	,3	,0	,2	,0	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,2	,1	,3	,1	,2	,2	,3	,2	,6	,5	,2	,6	,9	,8	,11	,9	,4	,9	,13	,14	,20	,10	,22	,26	,12	,25	,22	,19	,14	,19	,11	,21	,20	,14	,18	,15	,14	,19	,11	,7	,12	,11	,12	,12	,12	,4	,9	,5	,3	,1	,2	,2	,1	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,1	,0	,0	,4	,1	,0	,0	,2	,3	,4	,3	,4	,6	,4	,8	,7	,2	,2	,3	,11	,15	,14	,18	,14	,12	,24	,16	,26	,27	,16	,16	,23	,8	,8	,11	,14	,16	,24	,18	,24	,17	,12	,7	,22	,12	,8	,7	,6	,8	,4	,8	,3	,3	,1	,2	,0	,0	,2	,0	,1	,0	,0,
			0	,0	,1	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,1	,0	,0	,0	,3	,2	,1	,2	,2	,2	,5	,4	,7	,9	,4	,4	,4	,6	,10	,9	,16	,9	,12	,13	,20	,24	,18	,19	,27	,25	,19	,25	,25	,21	,16	,23	,19	,25	,20	,13	,20	,18	,15	,14	,14	,12	,6	,9	,9	,9	,2	,4	,2	,1	,2	,0	,0	,1	,0	,0	,0	,1,
			0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,1	,3	,3	,3	,1	,3	,6	,4	,5	,5	,4	,7	,5	,16	,11	,9	,19	,11	,12	,18	,16	,17	,22	,30	,27	,27	,28	,25	,29	,20	,37	,23	,16	,19	,13	,16	,16	,30	,21	,20	,21	,15	,10	,18	,7	,13	,7	,6	,4	,3	,2	,1	,1	,0	,0	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,3	,1	,3	,4	,3	,3	,4	,4	,9	,10	,8	,8	,14	,9	,16	,6	,11	,11	,19	,21	,22	,29	,24	,21	,16	,27	,30	,31	,26	,36	,35	,26	,24	,20	,34	,23	,24	,19	,17	,19	,18	,17	,12	,6	,5	,4	,5	,7	,2	,2	,3	,2	,2	,2	,0	,0	,0	,1,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,1	,0	,2	,3	,2	,2	,2	,1	,5	,4	,2	,3	,8	,6	,3	,7	,14	,11	,10	,13	,12	,19	,18	,19	,23	,21	,39	,24	,33	,22	,25	,29	,29	,38	,30	,29	,17	,25	,15	,12	,30	,27	,23	,18	,15	,15	,16	,12	,7	,7	,9	,4	,4	,3	,2	,0	,2	,0	,0	,1	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,2	,0	,6	,3	,8	,3	,3	,5	,9	,9	,10	,11	,13	,15	,13	,15	,16	,14	,23	,14	,21	,18	,27	,23	,30	,30	,29	,21	,35	,22	,31	,34	,25	,20	,32	,20	,21	,19	,22	,18	,15	,16	,12	,6	,7	,11	,8	,4	,3	,2	,0	,0	,1	,0	,0	,0	,2,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,1	,1	,1	,1	,3	,2	,10	,6	,4	,6	,10	,12	,10	,11	,13	,13	,13	,20	,23	,17	,24	,22	,26	,18	,25	,40	,28	,29	,42	,37	,35	,33	,42	,39	,30	,23	,25	,25	,19	,22	,16	,19	,14	,15	,14	,6	,11	,5	,3	,1	,2	,1	,2	,1	,0	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,2	,2	,3	,1	,1	,4	,2	,10	,2	,3	,7	,14	,13	,14	,12	,16	,21	,23	,22	,15	,21	,22	,30	,26	,30	,32	,18	,40	,27	,37	,37	,30	,34	,46	,29	,32	,34	,22	,34	,24	,38	,28	,22	,17	,14	,8	,8	,9	,5	,7	,5	,3	,2	,0	,2	,2	,0	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,1	,1	,0	,1	,0	,0	,1	,3	,2	,5	,2	,4	,5	,6	,6	,13	,12	,19	,9	,7	,21	,19	,18	,22	,29	,18	,17	,17	,27	,29	,37	,29	,27	,30	,43	,42	,37	,38	,50	,41	,46	,26	,26	,29	,30	,26	,15	,22	,22	,17	,19	,11	,6	,6	,6	,2	,4	,0	,1	,2	,1	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,2	,3	,0	,6	,8	,10	,19	,19	,13	,24	,17	,21	,16	,23	,23	,27	,23	,26	,33	,26	,47	,49	,38	,52	,39	,40	,40	,43	,34	,35	,40	,35	,36	,25	,27	,25	,25	,23	,16	,15	,16	,9	,7	,7	,9	,8	,4	,4	,1	,2	,0	,0	,1,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,1	,0	,0	,0	,0	,1	,1	,2	,8	,3	,4	,8	,6	,11	,11	,12	,15	,18	,13	,22	,22	,25	,31	,26	,35	,28	,29	,27	,31	,45	,51	,48	,40	,44	,55	,54	,32	,43	,47	,52	,30	,30	,26	,29	,26	,15	,19	,13	,6	,7	,11	,3	,3	,5	,6	,5	,0	,2	,0	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,1	,0	,1	,0	,0	,2	,2	,1	,3	,0	,1	,5	,6	,9	,17	,13	,20	,11	,18	,15	,24	,22	,23	,30	,25	,25	,24	,38	,31	,35	,31	,51	,42	,47	,44	,55	,48	,48	,39	,44	,41	,40	,32	,17	,21	,28	,22	,20	,11	,21	,16	,13	,6	,10	,6	,2	,4	,0	,1	,0	,2,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,2	,1	,0	,2	,3	,2	,4	,6	,5	,7	,10	,13	,15	,13	,22	,24	,21	,16	,20	,28	,30	,30	,32	,47	,29	,34	,37	,37	,43	,53	,46	,49	,49	,38	,38	,51	,36	,61	,34	,22	,21	,23	,26	,18	,19	,20	,12	,15	,3	,7	,9	,2	,1	,0	,2	,1	,0,
			0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,1	,0	,2	,0	,0	,0	,1	,1	,2	,0	,0	,2	,3	,4	,2	,7	,13	,10	,15	,15	,17	,24	,27	,24	,23	,26	,28	,24	,26	,42	,30	,29	,40	,40	,53	,46	,44	,54	,42	,50	,69	,50	,49	,38	,64	,44	,39	,30	,31	,23	,28	,23	,26	,13	,10	,7	,7	,3	,5	,4	,4	,3	,0,
			0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,2	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,1	,2	,2	,2	,2	,6	,4	,7	,9	,13	,19	,8	,17	,25	,29	,33	,36	,35	,41	,45	,39	,34	,40	,42	,43	,43	,51	,50	,42	,40	,45	,62	,56	,71	,54	,52	,49	,45	,27	,21	,26	,24	,21	,18	,12	,13	,9	,7	,5	,3	,1	,0	,1	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,1	,0	,0	,1	,0	,0	,0	,5	,1	,7	,5	,4	,8	,9	,7	,15	,23	,18	,17	,21	,26	,29	,28	,38	,43	,38	,42	,57	,42	,35	,45	,50	,52	,61	,51	,66	,54	,57	,43	,50	,61	,54	,40	,38	,22	,20	,20	,22	,25	,18	,10	,9	,5	,4	,2	,1	,2	,0	,1,
			0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,2	,2	,0	,1	,2	,1	,6	,3	,3	,4	,5	,11	,10	,20	,13	,25	,25	,32	,30	,35	,30	,40	,36	,43	,41	,55	,41	,47	,50	,54	,51	,70	,49	,57	,57	,62	,63	,72	,57	,59	,58	,55	,52	,24	,18	,22	,24	,21	,20	,7	,8	,6	,3	,4	,1	,2	,0,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,2	,0	,0	,0	,0	,1	,1	,1	,3	,2	,2	,9	,9	,13	,16	,21	,27	,29	,35	,34	,28	,35	,39	,37	,45	,74	,52	,41	,57	,48	,61	,51	,79	,65	,77	,62	,73	,60	,51	,45	,49	,48	,34	,39	,20	,17	,15	,19	,8	,12	,12	,8	,4	,4	,3	,0	,4,
			0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,2	,0	,0	,3	,3	,6	,4	,7	,8	,10	,13	,15	,30	,24	,26	,24	,46	,40	,37	,49	,47	,32	,53	,41	,65	,59	,58	,48	,66	,53	,53	,53	,68	,55	,64	,67	,60	,53	,47	,41	,42	,22	,25	,17	,21	,12	,11	,6	,4	,3	,4	,0	,2,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,1	,0	,0	,1	,0	,2	,4	,3	,3	,2	,6	,5	,2	,7	,6	,14	,7	,17	,26	,26	,25	,40	,30	,49	,41	,47	,59	,56	,49	,50	,64	,68	,64	,54	,74	,62	,69	,61	,59	,56	,72	,67	,53	,46	,37	,52	,39	,22	,14	,12	,13	,13	,10	,5	,7	,4	,1	,4,
			0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,1	,0	,0	,2	,0	,0	,2	,1	,1	,3	,1	,2	,3	,2	,3	,5	,7	,11	,13	,14	,13	,19	,18	,37	,32	,37	,38	,38	,55	,49	,57	,48	,54	,74	,62	,57	,57	,68	,61	,73	,60	,66	,69	,54	,61	,74	,61	,53	,42	,44	,33	,32	,21	,22	,13	,6	,10	,5	,9	,0	,1	,1,
			0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,2	,0	,0	,1	,1	,1	,3	,4	,5	,3	,2	,9	,3	,11	,16	,4	,18	,28	,31	,33	,31	,31	,40	,58	,60	,52	,59	,63	,62	,78	,52	,57	,67	,59	,75	,81	,69	,71	,64	,56	,56	,55	,47	,39	,41	,29	,37	,19	,15	,8	,5	,10	,7	,2	,6	,5,
			0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,2	,1	,0	,1	,0	,0	,4	,0	,0	,3	,1	,6	,2	,3	,5	,4	,9	,8	,16	,14	,27	,52	,49	,46	,39	,46	,41	,65	,45	,62	,76	,80	,58	,61	,51	,75	,67	,62	,85	,84	,65	,70	,55	,60	,58	,45	,47	,38	,36	,31	,9	,11	,10	,6	,3	,5	,3	,9)
col.names	<- c("5-9","10-14","15-19","20-24",as.character(seq(25,94)),"95+")
cases	<- matrix(data=v.cases,nrow=41, ncol=75, byrow=TRUE,dimnames=list(NULL,col.names))
if(all.age.groups==FALSE)
	return(apc.data.list(
				response		=cases[,seq(5,69)]	,
				data.format		="PA"			,
				age1			=25				,
				per1			=1967			,
				unit			=1				))
if(all.age.groups==TRUE)
	return(apc.data.list(
				response		=cases			,
				data.format		="PA"			,
				label			="UK mesothelioma"))
}	#	data.asbestos						

###############################
#	MOTOR DATA
###############################
data.loss.VNJ	<- function()
#	BN, 24 Apr 2015 (6 Feb 2015)
#	A Chain-Ladder with A,C effects
#
#	Taken from tables 1,2 of
#	Verrall R, Nielsen JP, Jessen AH (2010)
#	Prediction of RBNS and IBNR claims using claim amounts and claim counts
#	ASTIN Bulletin 40, 871-887
#
#	Also analysed in 
#
#	Martinez Miranda MD, Nielsen B, Nielsen JP and Verrall R (2011)
#	Cash flow simulation for a model of outstanding liabilities based on claim amounts and claim numbers
#	ASTIN Bulletin 41, 107-129
#
#	Kuang D, Nielsen B, Nielsen JP (2015)
#	The geometric chain-ladder
#	Scandinavian Acturial Journal, to appear
#
#	Data from Codan, Danish subsiduary of Royal & Sun Alliance 
#	Portfolio of motor policies: third party liability
#	Units in years
#	X Paid run-off triangle
#	N Number of reported claims
{	#	data.loss.VNJ
#	dimension
k		<- 10
#	Number of reported claims
Nvec	<- c(	6238  , 831   , 49    , 7     , 1 , 1 , 2 , 1 , 2 , 3 ,
				7773  , 1381  , 23    , 4     , 1 , 3 , 1 , 1 , 3 ,  
				10306 , 1093  , 17    , 5     , 2 , 0 , 2 , 2 ,    
				9639  , 995   , 17    , 6     , 1 , 5 , 4 ,    
				9511  , 1386  , 39    , 4     , 6 , 5 ,     
				10023 , 1342  , 31    , 16    , 9 ,  
				9834  , 1424  , 59    , 24    ,   
				10899 , 1503  , 84    ,       
				11954 , 1704  ,
				10989 )
#	X	as a vector
Xvec	<- c(   451288 , 339519 , 333371 , 144988 , 93243  , 45511  , 25217 , 20406 , 31482 , 1729 ,
			    448627 , 512882 , 168467 , 130674 , 56044  , 33397  , 56071 , 26522 , 14346 ,     
			    693574 , 497737 , 202272 , 120753 , 125046 , 37154  , 27608 , 17864 ,       
			    652043 , 546406 , 244474 , 200896 , 106802 , 106753 , 63688 ,       
			    566082 , 503970 , 217838 , 145181 , 165519 , 91313  ,      
			    606606 , 562543 , 227374 , 153551 , 132743 ,        
			    536976 , 472525 , 154205 , 150564 ,        
			    554833 , 590880 , 300964 ,          
			    537238 , 701111 ,           
			    684944 )     
return(c(apc.data.list(
			response	=vector.2.triangle(Xvec,k)	,
			data.format	="CL"						,
			time.adjust =1							,
			label		="loss VNJ"					),			
			list(
			counts		=vector.2.triangle(Nvec,k)	)))
}	#	data.loss.VNJ

###############################
#	LOSS TRIANGLE DATA
###############################
data.loss.BZ	<- function()
#	BN, 24 Apr 2014 (7 Feb 2015)
#	A Loss Triangle with A,P,C effects
#
#	Taken from table 3.5 of
#   Barnett, G. and Zehnwirth, B. (2000). Best estimates for reserves.
#	Proc. Casualty Actuar. Soc. 87, 245--321.
#
#	Also analysed in 
#
#	Kuang D, Nielsen B and Nielsen JP (2011)
#	Forecasting in an extended chain-ladder-type model
#	Journal of Risk and Insurance 78, 345-359
#
#	BZ write: "loss development array with a major trend change between payment years 1984 and 1985"
#	Time in years
#	X Paid run-off triangle
{	#	data.loss.BZ
#	dimension
k		<- 11
#	X	as a vector
Xvec	<- c(   153638,	188412,	134534,	 87456,	 60348,	42404,	31238,	21252,	16622,	14440,	12200,
			    178536,	226412,	158894,	104686,	 71448,	47990,	35576,	24818,	22662,	18000,        
			    210172,	259168,	188388,	123074,	 83380,	56086,	38496,	33768,	27400,                
			    211448,	253482,	183370,	131040,	 78994,	60232,	45568,	38000,                        
			    219810,	266304,	194650,	120098,	 87582,	62750,	51000,                                
			    205654,	252746,	177506,	129522,	 96786,	82400,                                        
			    197716,	255408,	194648,	142328,	105600,                                               
			    239784,	329242,	264802,	190400,                                                       
			    326304,	471744,	375400,                                                               
			    420778,	590400,                                                                       
				496200)
Exposure	<- c(     2.2,    2.4,	   2.2,	   2.0,	   1.9,	  1.6,    1.6,    1.8,    2.2,    2.5,    2.6)				
return(c(apc.data.list(
			response	=vector.2.triangle(Xvec,k)	,
			data.format	="CL"						,
			coh1		=1977						,
			time.adjust =1							,
			label		="loss BZ"					),
			exposure	=Exposure					))
}	#	data.loss.BZ

###############################
#	LOSS TRIANGLE DATA
###############################
data.loss.TA	<- function()
#	BN, 24 Apr 2015 (18 Mar 2015)
#	A Loss Triangle with A,C effects and over-dispersion
#
#	Attributed to
#	Taylor and Ashe
#
##	Analysed in
#
#	Verrall, R.J. (1991)
#	On the estimation of reserves from loglinear models
#	Insurance: Mathematics and Economics 10, 75-80
#
#	England, P., Verrall, R.J. (1999)
#	Analytic and bootstrap estimates of prediction errors in claims reserving
#	Insurance: Mathematics and Economics 25, 281-293
#
#	X Paid run-off triangle
{	#	data.loss.BZ
#	dimension
k		<- 10
#	X	as a vector
Xvec	<- c(	357848,	766940, 610542, 482940, 527326, 574398, 146342, 139950, 227229,  67948,
				352118, 884021, 933894,1183289, 445745, 320996, 527804, 266172, 425046,
				290507,1001799, 926219,1016654, 750816, 146923, 495992, 280405,
				310608,1108250, 776189,1562400, 272482, 352053, 206286,
				443160, 693190, 991983, 769488, 504851, 470639,
				396132, 937085, 847498, 805037, 705960,
				440832, 847631,1131398,1063269,
				359480,1061648,1443370,
				376686, 986608,
				344014)
return(apc.data.list(
			response	=vector.2.triangle(Xvec,k)	,
			data.format	="CL"						,
			time.adjust =1							,
			label		="loss TA"					))
}	#	data.loss.TA

###############################
#	AIDS reports in England and Wales
###############################
data.aids	<- function(all.age.groups=FALSE)
#	BN, 7 Feb 2016
#	Numbers of AIDS reports in England and Wales to the end of 1992 by quarter
#	
#	Attributed to
#	De Angelis and Gilks (1994)
#
#	Analysed in
#
#	Davison, A.C. and Hinkley, D.V. (1997) Bootstrap methods and their applications, Cambridge UP
{	#	data.aids
	#	data for coh 1983:3 to 1989:1
	v.cases		<- c(	2,	6,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,
						2,	7,	1,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
						4,	4,	0,	1,	0,	2,	0,	0,	0,	0,	2,	1,	0,	0,	0,
						0, 10,	0,	1,	1,	0,	0,	0,	1,	1,	1,	0,	0,	0,	0,
						6, 17,	3,	1,	1,	0,	0, 	0,	0,	0,	0,	1,	0,	0,	1,
						5, 22, 	1,	5,	2,	1,	0,	2,	1,	0,	0,	0,	0,	0,	0,
						4, 23, 	4,	5, 	2,	1,	3,	0,	1,	2,	0,	0,	0,	0,	2,
					   11, 11, 	6,	1,	1,	5,	0,	1,	1,	1,	1,	0,	0,	0,	1,
					    9, 22,	6,	2,	4,	3,	3,	4,	7,	1,	2,	0,	0,	0,	0,
						2, 28,	8,	8,	5,	2,	2,	4,	3,	0,	1,	1,	0,	0,	1,
						5, 26, 14,	6,	9,	2,	5,	5,	5,	1,	2,	0,	0,	0,	2,
						7, 49, 17, 11, 	4, 	7, 	5, 	7,	3,	1,	2,	2, 	0,	1,	4,
					   13, 37, 21,	9,	3,	5, 	7,	3,	1,	3,	1,	0,	0,	0,	6,
					   12, 53, 16, 21,	2,	7,	0,	7,	0,	0,	0,	0,	0,	1,	1,
					   21, 44, 29, 11,	6,	4,	2,	2,	1,	0,	2,	0,	2,	2, 	8,
					   17, 74, 13, 13,	3,	5,	3,	1,	2,	2,	0,	0,	0,	3,	5,
					   36, 58, 23, 14,  7,	4,	1,	2,	1,	3,	0,	0,	0,	3,	1,
					   28, 74, 23, 11,	8,	3,	3,	6,	2,	5,	4,	1,	1,	1,	3,
					   31, 80, 16,	9,	3,	2,	8,	3,	1,	4,	6,	2,	1,	2,	6,
					   26, 99, 27, 	9,	8, 11,	3,	4,	6,	3,	5,	5,	1,	1, 	3,
					   31, 95, 35, 13, 18,	4,	6,	4,	4,	3,	3,	2,	0,	3,	3,
					   36, 77, 20, 26, 11,	3,	8,	4,	8,	7,	1,	0,	0,	2,	2,
					   32, 92, 32, 10, 12, 19, 12,	4,	3,	2, 	0,	2,	2, 	0,	2,
					   15, 92, 14, 27, 22, 21, 12,	5,	3,	0,	3,	3,	0,	1,	1,
					   34,104, 29, 31, 18,	8,	6,	7, 	3,	8,	0,	2,	1,	2, NA, 
					   38,101, 34, 18, 	9, 15, 	6,	1,	2,	2,	2,	3,	2, NA, NA, 
					   31,124, 47, 24, 11, 15,	8,	6,	5,	3,	3,	4, NA, NA, NA, 
					   32,132, 36, 10,	9,	7,	6,	4,	4,	5,	0, NA, NA, NA, NA, 
					   49,107, 51, 17, 15,  8,	9,	2,	1,	1, NA, NA, NA, NA, NA, 
					   44,153, 41, 16, 11,	6,	5,	7,	2, NA, NA, NA, NA, NA, NA, 
					   41,137, 29, 33,	7, 11,	6,	4,	3, NA, NA, NA, NA, NA, NA, 
					   56,124, 39, 14, 12, 	7, 10,	1, NA, NA, NA, NA, NA, NA, NA, 
					   53,175, 35, 17, 13, 11, 	2, NA, NA, NA, NA, NA, NA, NA, NA, 
					   63,135, 24, 23, 12,	1, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
					   71,161, 48, 25, 	5, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
					   95,178, 39,	6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
					   76,181, 16, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
					   67, 66, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
col.names	<- c("0*",as.character(seq(1,13)),"14+")
cases.all	<- matrix(data=v.cases,nrow=38, ncol=15, byrow=TRUE,dimnames=list(NULL,col.names))
cases.clean	<- cases.all
for(row in 0:7)
	cases.clean[38-row,2+row]	<- NA
if(all.age.groups==FALSE)
	return(apc.data.list(
				response		=t(cases.clean)		,
				data.format		="trap"				,
				age1			=0					,
				coh1			=1983.5				,
				unit			=1/4				,
				per.zero		=0					,
				per.max			=38					,
				label			="UK AIDS - clean"	,
				))
if(all.age.groups==TRUE)
	return(apc.data.list(
				response		=cases.all			,
				data.format		="CA"				,
				age1			=0					,
				coh1			=1983.5				,
				unit			=1/4				,
				label			="UK AIDS - all: last column reporting delay >= 14, last diagonal: incomplete count",
				))	
}	#	data.aids
# 	apc.fit.table(data.Aids(),"od.poisson.response")
# 	fit <- apc.fit.model(data.Aids(),"poisson.response","AC")
#	forecast <- apc.forecast.ac(fit)
#	data.sums.coh <- apc.data.sums(data.Aids())$sums.coh
#	forecast.total <- forecast$response.forecast.coh
#	forecast.total[,1]	<- forecast.total[,1]+data.sums.coh[26:38]
#	plot(seq(1983.5,1992.75,by=1/4),data.sums.coh,xlim=c(1988,1993),ylim=c(200,600),main="Davison, Hinkley, Fig 7.6, parametric version")
#	apc.polygon(forecast.total,x.origin=1989.5,unit=1/4)