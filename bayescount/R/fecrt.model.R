#' @title Create an MCMC model to analyse FECRT Data
#' @name fecrt.model
#' @aliases fecrt.model FECRT.model
#' @export

#' @description
#' This function generates and compiles a JAGS model representation of a faecal egg count reduction test (FECRT) analysis.  The return value can be updated manually using \code{\link[runjags]{extend.jags}}

#' @details
#' Pre-treatment data are assumed to arise from either a gamma-Poisson or zero-inflated gamma-Poisson distribution, with post-treatment data described by a separate gamma-Poisson or zero-inflated gamma-Poisson distribution.  The change in mean between these distributions is therefore the mean egg count reduction.  A change in shape parameter of the gamma distribution is also permitted if fix.variation=FALSE. If paired.model=TRUE, a slightly different formulation is used whereby the observed count pre-treatment is assumed to follow a compound gamma-gamma-Poisson distribution with the variability within and between animals separated.  The post treatment mean for each animal is derived from the pre-treatment animal mean and FEC reduction.  This formulation allows data with non-random missing post-treatment counts to be analysed correctly, and also allows data with repeat counts from an individual to be analysed - providing a method of increasing the power of the method substantially.  The fix.efficacy=FALSE option is only permitted for the paired.model=TRUE option, and allows the reduction estimated to vary between individuals.  Problems with convergence are likely to be encountered with this option unless there is a substantial amount of replicate data avaialble within individuals.
#'

#' @keywords models

#' @return
#' Returns an object of class \code{\link[runjags]{runjags-class}}

#' @seealso
#' \code{\link{fecrt.analysis}} for comparisons of fitted MCMC models to bootstrapping results

#' @references
#' M. J. Denwood, S. W. J. Reid, S. Love, M. K. Nielsen, L. Matthews, I. J. McKendrick, and G. T. Innocent. Comparison of three alternative methods for analysis of equine Faecal Egg Count Reduction Test data. Prev. Vet. Med. (2009), doi:10.1016/j.prevetmed.2009.11.009

#' @examples
#' \dontrun{
#' # Data in an appropriate format:
#' data <- data.frame(Count=rpois(80,rep(c(10,10,2,2), 20)), 
#' Subject=rep(1:20, each=4), Time=rep(rep(1:2,each=2),40), 
#' Sample=1:2, Control=rep(c(0,1), each=40))
#' # Compile the model - a paired model is required because 
#' # there are replicate samples within an individual:
#' model <- fecrt.model(data, paired.model=TRUE)
#' # Update the model - requires runjags:
#' library('runjags')
#' results <- extend.jags(model, burnin=5000)	
#' }	

#' @param data a data frame containing the dataset to be analysed in long format, with columns of Count (the number of eggs observed), Subject (a name or number uniquely identifiyig each animal/patient) and Time (1 = Pre-treatment, 2 = Post-treatment), and optional columns of Sample (a number) and Control (1 = Treatment group, 2 = Control group).  Multiple Counts are permitted from the same subject, either with the same Sample number (repeated count observations from the same processed laboratory sample) or different Sample numbers (independently analysed samples from the same individual).

#' @param paired.model the FECRT can be run using two compound gamma distributions to describe the variability within and between animals in place of the single gamma distribution combining the two sources of variability.  When using two compound gamma distributions for pre-treatment data, the post-treatment data are paired to the pre-treatment data by animal.  The advantage of using a single distribution is that the model is more identifiable, and therefore is likely to converge more quickly and return errors less frequently.  The advantage of the paired model is that it allows repeat measurements tobe incorporated in the model, and non-randomly missing data (ie.protocols that involve post-treatment sampling of only animals with a high pre-treatment count) to be modelled appropriately.  The simple model cannot be used when using repeat samples within animal, and will provide inaccurate results if animals are targeted for post-treatment sampling based on their pre-treatment count.  Default uses the simple model (FALSE).

#' @param fix.controls the FECRT can be run either using control animals to estimate the overall reduction atrributable to the treatment, or assuming that the true mean of the controls does not change over time and only using this data to improve the estimate of the pre-treatment mean.  The former is more in the spirit of what is intended with controls, but the latter will give smaller confidence intervals at the cost of the strong assumption that the underlying population mean has not changed between the treatment dates.  Default FALSE.

#' @param fix.efficacy force all treatment animals to have the same reduction following treatment, so that all residual variation is absorbed by the extra-Poisson sampling variation. Setting this to FALSE requires the paired model and a substantial number of replicate samples within subjects.  Default TRUE.

#' @param fix.variation force the post-treatment extra-Poisson variation to be the same as that pre-treatment.  This parameter is difficult to estimate, so forcing it to be the same as the pre-treatment value may improve convergence, although it is a strong assumption.  Note that this parameter and the efficacy variance parameter are likely to be strongly inter-dependent.  Default FALSE.

#' @param zero.inflation option to use a zero-inflated gamma Poisson inplace of the gamma Poisson distribution.  If TRUE, zero inflated distributions are used pre- and post-treatment (with the prevalence fixed between pre- and post-treatment distributions).  Default FALSE.

#' @param effect.prior the prior distribution to use for the change in mean egg count.  Any syntactically valid distribution in JAGS may be used, but it must include the default initial values of 0.01 and 1. Default is a uniform distribution between 0 and 1.

#' @param mean.prior the prior distribution to use for the pre-treatment mean egg count.  Any syntactically valid strictly positive distribution in JAGS may be used, but it must include the default initial value of the arithmetic mean pre-treatment count. Default is DuMouchel's prior distribution.

#' @param precision.prior the prior distribution to use for the shape parameter (k) of the various gamma distributions used.  Any syntactically valid strictly positive distribution in JAGS may be used, but it must include the default initial values of 0.1 and 10. Default is DuMouchel's prior distribution.

#' @param ... other options to be passed directly to \code{\link[runjags]{autorun.jags}} when setting up the model.  This may include initial values, for example .RNG.seed and .RNG.name to allow the analysis to be reproducible.


#' @rdname fecrt.model
fecrt.model <- function(data, paired.model=FALSE, fix.controls=FALSE, fix.efficacy=TRUE, fix.variation=TRUE, zero.inflation=FALSE, effect.prior='dbeta(1,1)', mean.prior='dmouch(1)', precision.prior='dmouch(1)', ...){
	# check runjags version 2
	if(!utils::packageDescription('runjags', field='Version') >= 2)
		stop('The fecrt.model function requires version 2 of the runjags package', call.=FALSE)
	
	data <- sortfecrtdata(data)
	jagsdata <- as.list(data)
	jagsdata$SubjectNames <- NULL
	
	if(all(jagsdata$Control==1)){
		fix.controls <- TRUE
	}
	
  	if(any(jagsdata$Sample>1) && !paired.model)
		stop('Replicate samples are only possible with paired.model=TRUE', call.=FALSE)
	if(!fix.efficacy && !paired.model)
		stop('The fix.efficacy=FALSE option is only available with paired.model=TRUE', call.=FALSE)
	
	names <- data$SubjectNames
	
	passed <- list(...)

	if(any(names(passed)=='n.chains')){
		n.chains <- passed$n.chains
		if(n.chains!=2)
			stop('The fecrt.model requires exactly 2 chains', call.=FALSE)		
	}else{
		n.chains <- 2
	}
	
	if(any(names(passed)=='inits')){
		inits <- passed$inits
		if(!is.list(inits))
			stop('The initial values must be provided as a list of named lists')
		if(!is.list(inits[[1]]))
			inits <- lapply(1:n.chains, function(x) return(inits))
		if(!is.na(n.chains) && length(inits)!=n.chains)
			stop('Invalid initial variables provided - this must be a list of length matching the n.chains provided', call.=FALSE)
	}else{
		inits <- lapply(1:n.chains, function(x) return(list()))
	}
	
	# not allowed in ...:  data, mutate
	notallowed <- c('data','mutate','monitor')
	if(any(notallowed%in%names(passed)))
		warning(paste('The following arguments are not allowed for fecrt.model and have been : ', paste(notallowed, collapse=', '), sep=''))
	passed <- passed[!names(passed)%in%c('inits','n.chains',notallowed)]
	
	replicates <- FALSE
	NCounts <- length(jagsdata$Count)
	NSubjects <- max(jagsdata$Subject)
	NPre=NPost <- numeric(NSubjects)
	for(i in 1:NSubjects){
		NPre[i] <- max(data$Sample[data$Subject==i & data$Time==1], 0)
		NPost[i] <- max(data$Sample[data$Subject==i & data$Time==2], 0)
	}
	stopifnot(all(!is.na(NPre)) && all(!is.na(NPost)))
	jagsdata <- c(jagsdata, list(NCounts=NCounts, NSubjects=NSubjects, NPre=NPre, NPost=NPost, MaxPre=max(NPre), MaxPost=max(NPost)))
	
	model <- paste("
model{
	for(i in 1:NCounts){
		Count[i] ~ dpois(lambda[Subject[i],Time[i],Sample[i]])
	}
	for(a in 1:NSubjects){
		for(s in 1:NPre[a]){
			lambda[a,1,s] <- animalmean[a,1]", if(zero.inflation) " * non_zero[a]", " * presamplegamma[a,s]
		}
		for(s in 1:NPost[a]){
			lambda[a,2,s] <- animalmean[a,2] ", if(zero.inflation) " * non_zero[a]", " * postsamplegamma[a,s]
		}
		for(s in 1:MaxPre){
			presamplegamma[a,s] ~ dgamma(presample_k, presample_k)T(10^-200,)
		}		
		for(s in 1:MaxPost){
			postsamplegamma[a,s] ~ dgamma(postsample_k, postsample_k)T(10^-200,)
		}		
		animalmean[a,1] <- groupmean ", if(paired.model) paste("* animalgamma[a]
		animalgamma[a] ~ dgamma(animal_k, animal_k)T(10^-200,)
		animalmean[a,2] <- animalmean[a,1] * deltamean[a,Control[a]]
		deltamean[a,2] <- meanchange[2]
		deltamean[a,1] <- meanchange[1]", if(!fix.efficacy) " * efficacygamma[a]
		efficacygamma[a] ~ dgamma(efficacy_k, efficacy_k)T(10^-200,)", "
		", sep="") else paste("
		animalmean[a,2] <- groupmean * meanchange[Control[a]]", sep=""), "
		", if(zero.inflation) "non_zero[a] ~ dbern(nonzeroprob)", "
	}
	
	groupmean ~ ", mean.prior, "
	meanchange[1] ~ ", effect.prior, "
	", if(!fix.controls) paste("meanchange[2] ~ ", effect.prior, sep="") else "meanchange[2] <- 1", "
	", if(zero.inflation) "\nnonzeroprob ~ dbeta(1,1)", "
	", if(paired.model) paste("animal_k ~ ", precision.prior, sep=""), "
	presample_k ~ ", precision.prior, "
	postsample_k", if(fix.variation) " <- presample_k" else paste(" ~ ", precision.prior, sep=""), "
	", if(!fix.efficacy) paste("efficacy_k ~ ", precision.prior, "\n", sep=""), "	
	#monitor# groupmean, meanchange[1], ", if(!fix.controls) "meanchange[2], ", "presample_k, ", if(paired.model) "animal_k, ", if(!fix.variation) "postsample_k, ", if(!fix.efficacy) "efficacy_k, ", if(zero.inflation) "nonzeroprob, ", "
}		
", sep="")
# Dont specify the runjags module - it will be picked up automatically by runjags unless paretoprior is specified
# use glm module???  Doesn't seem to do any harm, but doesn't make a difference - all updated sequantially
	
	stopifnot(all(c('Count','Sample','Subject','Time')%in%names(jagsdata)))
	stopifnot(any(NPre>0) && any(NPost>0))
	stopifnot(NCounts > 0)
	
	newinits <- sapply(1:2, geninits, jagsdata=jagsdata, paired.model=paired.model, fix.controls=fix.controls, fix.efficacy=fix.efficacy, fix.variation=fix.variation, zero.inflation=zero.inflation, passedinits=inits)
	
	args <- passed
	args$model <- model
	args$data <- dump.format(jagsdata)
	args$inits <- newinits
	args$n.chains <- n.chains
	args$startsample <- 0
	args$startburnin <- 0
	args$adapt <- 0
	
	if(fix.controls){
		args$mutate <- function(x){
			reduction <- 100 * (1 - x[,'meanchange[1]',drop=FALSE])
			dimnames(reduction) <- list(dimnames(reduction)[[1]], 'reduction(%)')
			return(coda::as.mcmc(reduction))
		}
	}else{
		args$mutate <- function(x){
			reduction <- 100 * (1 - x[,'meanchange[1]',drop=FALSE]/x[,'meanchange[2]',drop=FALSE])
			dimnames(reduction) <- list(dimnames(reduction)[[1]], 'reduction(%)')
			return(coda::as.mcmc(reduction))
		}
	}
	
	# This is autorun.jags so it matches arguments that may be passed to fecrt.analysis:
	rjo <- do.call('autorun.jags', args)
	
	cat('Returning compiled model\n')
	return(rjo)

}



sortfecrtdata <- function(data){
	
	if('Animal'%in%names(data) && !'Subject'%in%names(data))
		data$Subject <- data$Animal
	
	if(!is.data.frame(data) || !all(c('Count','Subject', 'Time') %in% names(data)))
		stop('Invalid data format - it must be provided as a data frame of Count, Subject and Time, with optional Sample and Control - consult the help file for ?fecrt.model', call.=FALSE)
	if(is.null(data$Sample))
		data$Sample <- 1
	if(is.null(data$Control))
		data$Control <- 0   # Will add 1 below
	
	if(all(data$Control==2))
		stop('No treatment subjects found', call.=FALSE)
	if(all(data$Control==1))
		data$Control <- 0   # Will add 1 below
	
	if(any(is.na(data$Sample) | is.na(data$Subject) | is.na(data$Time)))
		stop('Missing values not allowed in Sample, Subject or Time', call.=FALSE)
	
	data$Control <- as.numeric(data$Control)
	if(all(data$Control %in% c(0,1)))
		data$Control <- data$Control + 1
	if(!all(data$Control %in% c(1,2)))
		stop('The Control element of the data frame must either be a numeric/logical of 0 (treatment) and 1 (control), or a numeric or factor with levels 1 (treatment) and 2 (control)', call.=FALSE)

	data$Time <- as.numeric(data$Time)
	if(all(data$Time %in% c(0,1)))
		data$Time <- data$Time + 1
	if(!all(data$Time %in% c(1,2)))
		stop('The Time element of the data frame must either be a numeric/logical of 0 (pre-treatment) and 1 (post-treatment), or a numeric or factor with levels 1 (pre-treatment) and 2 (post-treatment)', call.=FALSE)
	if(!any(data$Time==1))
		stop('No pre-treatment observations', call.=FALSE)
	if(!any(data$Time==2))
		stop('No post-treatment observations', call.=FALSE)
	
	if(any(is.na(data$Control[data$Time==2])))
		stop('Missing Control information is not allowed for the post-treatment data (Control is ignored for pre-treatment observations)', call.=FALSE)
	data$Control[data$Time==1] <- NA
	
	# Control is animal-centric so shorter
	data <- as.list(data)
	nc <- sapply(unique(data$Subject), function(a){
		ctls <- as.numeric(na.omit(data$Control[data$Subject==a]))
		if(length(ctls)==0)
			return(1)  # for animals with missing post-tx data it doesn't matter
		if(!all(ctls==ctls[1]))
			stop(paste('Inconsistent Control variable for subject ', a, sep=''))
		return(ctls[1])		
	})
	data$Control <- nc
	data$Subject <- factor(data$Subject)
	data$SubjectNames <- levels(data$Subject)
	data$Subject <- as.numeric(data$Subject)
	
	data <- data[names(data)%in%c('Count','SubjectNames', 'Subject', 'Time' ,'Sample','Control')]
	return(data)
}


geninits <- function(chain, jagsdata, paired.model, fix.controls, fix.efficacy, fix.variation, zero.inflation, passedinits=list()){
	
	if(class(passedinits)=='list')
		passedinits <- passedinits[[chain]]
	if(class(passedinits)%in%c('character','runjagsdata','runjagsinits'))
		passedinits <- list.format(passedinits[chain])
	
	premeans <- matrix(NA, nrow=max(jagsdata$Subject), ncol=max(jagsdata$NPre))
	postmeans <- matrix(NA, nrow=max(jagsdata$Subject), ncol=max(jagsdata$NPost))

	for(i in 1:ncol(premeans)) premeans[,i] <- sapply(1:nrow(premeans), function(x) return(if(!any(jagsdata$Subject==x & jagsdata$Time==1 & jagsdata$Sample==i)) NA else mean(jagsdata$Count[jagsdata$Subject==x & jagsdata$Time==1 & jagsdata$Sample==i])))

	for(i in 1:ncol(postmeans)) postmeans[,i] <- sapply(1:nrow(postmeans), function(x) return(if(!any(jagsdata$Subject==x & jagsdata$Time==2 & jagsdata$Sample==i)) NA else mean(jagsdata$Count[jagsdata$Subject==x & jagsdata$Time==2 & jagsdata$Sample==i])))

	animalpremeans <- sapply(unique(jagsdata$Subject), function(x) return(if(!any(jagsdata$Subject==x & jagsdata$Time==1)) NA else mean(jagsdata$Count[jagsdata$Subject==x & jagsdata$Time==1])))
	
	groupmean <- mean(premeans,na.rm=TRUE)

	stopifnot(chain%in%1:2)
	if(chain==1){
		# Chain 1:  low zi, high dispersion, high efficacy
		meanchange <- c(0.01, 1)
		presample_k <- 0.1
		postsample_k <- 0.1
		animal_k <- 0.1
		efficacy_k <- 0.1
		nonzeroprob <- 0.05
		non_zero <- rep(1, max(jagsdata$Subject)) 
	}else{
		# Chain 1:  high zi, low dispersion, low efficacy
		meanchange <- c(0.5, 1)
		presample_k <- 10
		postsample_k <- 10
		animal_k <- 10
		efficacy_k <- 10
		nonzeroprob <- 0.95		
		# It is fine for postmeans to be all missing - just evals as FALSE:
		non_zero <- as.numeric(apply(premeans>0,1,any,na.rm=TRUE) | apply(postmeans>0,1,any,na.rm=TRUE))
	}

	if(paired.model){
		deltamean <- pmax((apply(postmeans,1,sum,na.rm=TRUE) / apply(premeans,1,sum,na.rm=TRUE)), 0.001, na.rm=TRUE)
		deltamean[deltamean==Inf] <- 10
		# Missing will be NA which is fine:
		animalgamma <- pmax(animalpremeans/groupmean,0.1)
	}else{
		deltamean <- mean(postmeans,na.rm=TRUE) / mean(premeans,na.rm=TRUE)
		if(deltamean==Inf)
			deltamean <- 10
		animalgamma <- 1
	}
	deltamean <- pmax(deltamean, 0.01)
	
	efficacygamma <- deltamean / meanchange[1]
	
	# Missing will be NA which is fine:
	presamplegamma <- pmax(premeans/(groupmean*animalgamma),0.01)
	# Missing will be NA which is fine:
	postsamplegamma <- pmax(postmeans/(groupmean*animalgamma*deltamean),0.1)
	
	if(fix.controls)
		meanchange[2] <- NA

	toret <- list(groupmean=groupmean, meanchange=meanchange, presample_k=presample_k, presamplegamma=presamplegamma, postsamplegamma=postsamplegamma)
	if(paired.model)
		toret <- c(toret, list(animalgamma=animalgamma, animal_k=animal_k))
	if(!fix.efficacy)  # only possible with the paired model
		toret <- c(toret, list(efficacygamma=efficacygamma, efficacy_k=efficacy_k))
	if(!fix.variation)   	# only possible with the paired model
		toret <- c(toret, list(postsample_k=postsample_k))
	if(zero.inflation)
		toret <- c(toret, list(nonzeroprob=nonzeroprob, non_zero=non_zero))
	
	# Overwrite with anything we are given:
	for(n in names(passedinits))
		toret[[n]] <- passedinits[[n]]
	
	inits <- dump.format(toret)
	class(inits) <- 'runjagsinits'
	return(inits)
}

FECRT.model <- fecrt.model
