# A FUNCTION RETURNING THE LIKELIHOOD VALUE
# SO USER CAN USE THEIR OWN OPTIMISER
NB.likelihood<-function(N, infile, alleles, sample.interval)
{
# INPUT CHECKING
if (length(sample.interval)<2)
	{stop('sample.interval must be a vector of length 2 or more. See help for details.')}
if (any(sample.interval<0))
	{stop('sample.interval must be a non-negative vector. See help for details.')}
if (any(sample.interval%%1!=0))
	{stop('sample.interval must contain non-negative integers. See help for details.')}
if (any(alleles%%1!=0))
	{stop('alleles must be positive integers. See help for details.')}

	#####
	# INFILE TOOL
	dirmulti.infile<-function(infile, alleles, sample.interval)
	{
		f<-function(dat.line)
		{
		temp<-strsplit(dat.line, ' ')[[1]]
		temp<-as.numeric(temp[temp!=''])
		return(temp)
		}
	
	dat<-readLines(infile)
	dat<-sapply(dat, f)
	dat.length<-sapply(dat, length)
	dat<-dat[dat.length!=0]
	id<-(1:length(dat))%%length(alleles)
	
	temporal.samples<-length(sample.interval)
	
	# INFILE DIMENSION CHECKING
	if (length(dat)!=temporal.samples*length(alleles)) 
		{stop('infile format. See help for details.')}
	
	output<-list()
	
	for (i in 1:length(alleles))
		{output[[i]]<-matrix(unlist(dat[id==(i-1)]), ncol=alleles[i], byrow=T)}
	return(output)
	}
	#####
	# DIRICHLET-MULTINOMIAL PMF
	ddirmulti<-function(x, alpha, log=T)
	{
	if (any(x%%1!=0))
		{stop('allele counts must be non-negative integers')}
	if (any(x<0))
		{stop('allele counts must be non-negative integers')}
	temp<-lgamma(sum(x)+1)-sum(lgamma(x+1))+lgamma(sum(alpha))-lgamma(sum(alpha)+sum(x))+sum(lgamma(alpha+x)-lgamma(alpha))
	if (log==T) 
		{return(temp)}
	else 
		{return(exp(temp))}
	}
	#####
	# PARAMETER UPDATING TOOL
	dirichlet.updating<-function(dat.list, N.dip, sample.interval=sample.interval)
	{
	dirichlet.parms<-dat.list*0
	dirichlet.parms[1,]<-1
	for (i in 2:temporal.samples)
		{
		time.diff<-sample.interval[i]-sample.interval[i-1]
		kt<-(1-1/(2*N.dip))^time.diff
		drift.parms<-kt/(1-kt)
		temp<-dirichlet.parms[i-1,]+dat.list[i-1,]
		dirichlet.parms[i,]<-temp*drift.parms/(1+sum(temp)+drift.parms)
		}
	return(dirichlet.parms)
	}
	#####
	# LIKELIHOOD FOR EACH LOCUS
	dirichlet.log.likelihood<-function(dat.list, N.dip, sample.interval)
	{
	dirichlet.parms<-dirichlet.updating(dat.list=dat.list, N.dip=N.dip, sample.interval=sample.interval)
	likelihood.value<-rep(0, nrow(dat.list))
	for (i in 2:temporal.samples) 
		{likelihood.value[i]<-ddirmulti(x=dat.list[i,], alpha=dirichlet.parms[i,])}
	return(sum(likelihood.value))
	}
	#####
	# NEED TO WRAP IT ACROSS LOCI
	lapply.wrapper<-function(N.dip, dat, z=0)
	{
	log.likelihood.overall<-lapply(dat, dirichlet.log.likelihood, N.dip=N.dip, sample.interval)
	return(sum(unlist(log.likelihood.overall))-z)
	}
dat<-dirmulti.infile(infile=infile, sample.interval=sample.interval, alleles=alleles)
temporal.samples<-length(sample.interval)
log.likelihood<-lapply.wrapper(N.dip=N, dat=dat)
return(log.likelihood)
}
