setClass("homophily.estimate",
		representation(estimate="table",
				outcome.variable="character",
				weight.type="character",
				sample.mean="numeric",
				recruitment="logical",
				mean.group.weights="numeric"))


setMethod("initialize",
		"homophily.estimate",
		function(.Object,
				estimate,
				outcome.variable,
				weight.type,
				recruitment) {
			
			
			.Object@estimate <- estimate
			.Object@outcome.variable <- outcome.variable
			.Object@weight.type <- weight.type
			.Object@recruitment <- recruitment
			
			return(.Object)
		})

# A "friendly" constructor.  

homophily.estimate <- function(estimate,
		outcome.variable,
		weight.type,
		recruitment){
	new("homophily.estimate", estimate,
			outcome.variable,
			weight.type,
			recruitment)
}

setMethod("show",signature="homophily.estimate",
		definition=function(object){
			if(object@recruitment){
				cat("Recruitment Homophily for",object@outcome.variable,"\n\n")
				cat("Homophily =",attr(object@estimate,"homophily"),"\n\n")
				print(object@estimate)
				cat("\n")
				print(summary.h.table(object@estimate))
				cat("\n")
			}else{
				cat("Population Homophily Estimate for",object@outcome.variable,"\n")
				if(length(object@estimate)>1){
					print(object@estimate)
				}else{
					cat(format(object@estimate,nsmall=3),fill=TRUE)
				}
			}
		}
)


################################################################################
# A function for computing interval Estimates.  Do not export the local function.   
################################################################################


homophily.estimates.local <- function(rds.data,outcome.variable,
		weight.type,uncertainty,recruitment,N,to.group0.variable,to.group1.variable,subset,
		number.ss.samples.per.iteration, confidence.level){
	
	# The purpose of rds.data.frame is that it enforces the required
	# consistency on our data frames -- namely there must be a recruitment
	# structure described by 'id' and 'recruiter.id' fields, all non-seed
	# respondents must have a recruiter.   This means that if we have a 
	# valid rds.data.frame object, then we're free to convert it back to a 
	# simple (and slightly smaller) data frame for easy manipulation.
	
	
	if(is(rds.data,"rds.data.frame")){
		
		# Stop right now and cause an error if the group
		# variable does not appear in the data frame.
		
		if(!(outcome.variable %in% names(rds.data))){
			stop(sprintf("No variable called %s appears in the data.",
							outcome.variable))
		}
		
		if(!is.null(to.group0.variable)) stopifnot(to.group0.variable %in% names(rds.data))
		if(!is.null(to.group1.variable)) stopifnot(to.group1.variable %in% names(rds.data))
		if(!is.null(to.group0.variable)) {
			temp.data.frame <- data.frame(outcome.variable=as.vector(rds.data[[outcome.variable]]),
					to.group0.variable=as.vector(rds.data[[to.group0.variable]]),
					to.group1.variable=as.vector(rds.data[[to.group1.variable]]))
			colnames(temp.data.frame)[(1:length(outcome.variable))+(1:2)] <- 
					c(to.group0.variable, to.group1.variable)
			deg0 <- temp.data.frame[[to.group0.variable]]	
			deg1 <- temp.data.frame[[to.group1.variable]]	
		}else{
			temp.data.frame <- data.frame(outcome.variable=as.vector(rds.data[[outcome.variable]]))
		}
		colnames(temp.data.frame)[1:length(outcome.variable)] <- outcome.variable
		temp.data.frame$id <- as.character(rds.data[[attr(rds.data,"id")]])
		temp.data.frame$recruiter.id <- as.character(rds.data[[attr(rds.data,"recruiter.id")]])
		temp.data.frame$network.size <- rds.data[[attr(rds.data,"network.size.variable")]]
		temp.data.frame$wave <- rds.data$wave
		attr(temp.data.frame,"network.size.variable") <- "network.size"
		network.size <- "network.size"
		rds.data <- as.rds.data.frame(temp.data.frame,population.size=attr(rds.data,"population.size.mid"))
		
	}
	else{
		stop("The data is not an rds.data.frame object")
		
	}
	
	if(is.null(N)){
		N <- attr(rds.data, "population.size.mid")
	}
	if(is.null(weight.type)){
		if(is.null(N))
			weight.type <- "RDS-II"
		else
			weight.type <- "Gile's SS"
	}
	weight.type <- match.arg(weight.type,
			c("Gile's SS","RDS-I", "RDS-II", "RDS-I (DS)","Good-Fellows","Arithmetic Mean"))
	if(is.na(weight.type)) { # User typed an unrecognizable name
		stop(paste('You must specify a valid weight.type. The valid types are "Gile\'s SS","RDS-I", "RDS-II", "RDS-I (DS)", "Good-Fellows", and "Arithmetic Mean"'), call.=FALSE)
	}
	if(is.null(N)  && weight.type=="Gile's SS" && ! recruitment){
		stop("Parameter N missing, with no default for this data set")
	}
	
	if(is.null(number.ss.samples.per.iteration)){
		number.ss.samples.per.iteration <- 5000
	}
	
	
	if(is.null(confidence.level)){
		confidence.level <- 0.95
	}
	
	remvalues <- rds.data[[network.size]]==0 | is.na(rds.data[[network.size]])
	if(any(remvalues)){
		warning(paste(sum(remvalues),"of",nrow(rds.data),
						"network sizes were missing or zero. The estimator will presume these are",max(rds.data[[network.size]],na.rm=TRUE)), call. = FALSE)
		rds.data[[network.size]][remvalues] <- max(rds.data[[network.size]],na.rm=TRUE)
	}
	
	rds.data.nomiss <- rds.data
	rds.data.sub <- rds.data.nomiss
	deg <- rds.data.nomiss[[network.size]]
	
	outcome <- as.vector(rds.data.sub[[outcome.variable]])
	id <- as.character(rds.data.sub[["id"]])
	rid <- as.character(rds.data.sub[["recruiter.id"]])
	
	if(recruitment){
#	  So compute the recruitment homophily
		group.memberships <- as.vector(rds.data.nomiss[[outcome.variable]])
		
		group.names <- unique(group.memberships)
		
		# NB: if one wants to treat missing values as a separate group,
		# then the following will need to be changed. 
		group.names <- group.names[!is.na(group.names)]
		
		number.of.groups <- length(group.names)
		
		# Now translate these to indices.  This is done for the sake of efficiency.
		
		group.indices <- match(group.memberships,group.names)
		
		# We also want to extract the wave information.  We will do this in order to 
		# sort the observations by wave.
		wave <- get.wave(rds.data.nomiss)
		
		# If everything is ordered by ascending wave we can
		# simulate by running through the data set by rows.
		
		wave.order <- order(wave)    
		
		group.indices <- group.indices[wave.order]
		id <- as.character(rds.data.nomiss[["id"]])[wave.order]
		recruiter.id <- as.character(rds.data.nomiss[["recruiter.id"]])[wave.order]
		degrees <- data.frame(rds.data.nomiss)[wave.order,network.size]    
		wave <- wave[wave.order]
		
		# Now rationalize the recruiter.id information so that it corresponds to recruiter
		# row.
		
		recruiter.row <- match(recruiter.id,id)
		
		# The zeros in the recruiter id will be mapped to NA's
		recruiter.row[is.na(recruiter.row)] <- 0
		
		wave.one.start <- match(1,wave)
		seed.rows <- 1:(wave.one.start - 1)
		number.of.seeds <- length(seed.rows)
		sample.size <- length(id)
		
		
		# Now we need a count of transitions.
		tij <- matrix(0,nrow=number.of.groups,ncol=number.of.groups)
		
		for(row.index in wave.one.start:sample.size)
		{
			recruit.group.index <- group.indices[row.index]
			
			recruiter.group.index <- group.indices[recruiter.row[row.index]]
			
			if( (length(recruiter.group.index)>0) && (length(recruit.group.index)>0) && (!is.na(recruit.group.index)) && (!is.na(recruiter.group.index)))
			{
				tij[recruiter.group.index, recruit.group.index] <- 1 + tij[recruiter.group.index, recruit.group.index]
			}
			
		}
		
		colnames(tij) <- group.names
		rownames(tij) <- group.names
		names(dimnames(tij)) <- c(paste(outcome.variable," of respondent",sep=""),
				paste(outcome.variable," of recruit",sep=""))
		class(tij)<-"table" 
		
		homophily.group <- diag(tij)/ (apply(tij,1,sum)^2/sum(tij))
		num.grp <- rep(0,length=length(group.names))
		for(i in seq(along=group.names)){
			num.grp[i] <-	mean(group.memberships==group.names[i],na.rm=TRUE)
		}
# Original definition:
#		attr(tij,"homophily") <- sum(num.grp*homophily.group,na.rm=TRUE)
# Yuyan's new definition:
		s=sum(apply(tij,1,sum)*apply(tij,2,sum))/sum(tij)
		attr(tij,"homophily") <- sum(diag(tij))/s

		
		estimate <- tij 
	}else{
#	  So compute the population homophily
		weights.nomiss <- compute.weights(rds.data.sub,
				weight.type=weight.type,
				N=N,
				outcome.variable=outcome.variable,
				number.ss.samples.per.iteration=number.ss.samples.per.iteration,
				hajek=TRUE)
		
		weights.all <- rep(0,length=nrow(rds.data.nomiss))
		weights.all[subset] <- weights.nomiss
		outcome <- as.vector(rds.data.nomiss[[outcome.variable]])
		if(length(unique(outcome)) == 2){
			outcome <- as.numeric(outcome==max(outcome,na.rm=TRUE))
			
			weightsi <- weights.nomiss
			if(is.null(to.group0.variable)) {
				deg0 <- rep(NA,length(deg))
				deg1 <- rep(NA,length(deg))
				for(k in 1:nrow(rds.data.nomiss)){
					o=outcome[rid == id[k]]
					if(length(o)>0){
						deg0[k] = sum(!o) * deg[k] / length(o)
						deg1[k] = sum( o) * deg[k] / length(o)
					}else{
						weightsi[k] <- 0
					}
				}
				weightsi[is.na(deg0)|is.na(deg1)] <- 0
				deg0[is.na(deg0)] <- 0
				deg1[is.na(deg1)] <- 0
			}else{
				deg0 <- as.vector(rds.data.nomiss[,to.group0.variable])	
				deg1 <- as.vector(rds.data.nomiss[,to.group1.variable])
			}
			prop.outcome <- HT.estimate(weights=weightsi,outcome=outcome)
			mean.degree.outcome0 <- HT.estimate(weights=weightsi,outcome=deg*(outcome==0)) / (1-prop.outcome)
			mean.degree.outcome1 <- HT.estimate(weights=weightsi,outcome=deg*(outcome==1)) / prop.outcome
			mean.degree <- HT.estimate(weights=weightsi,outcome=deg)
			s01 <- (outcome==1)*deg0 + (outcome==0)*deg1
			mean.s01 <- HT.estimate(weights=weightsi,outcome=s01)
			
			estimate <- 2*prop.outcome*mean.degree.outcome0*(1-prop.outcome)*mean.degree.outcome1 / (mean.s01*mean.degree)
		}else{
			weights.nomiss <- compute.weights(rds.data,
					weight.type=weight.type, N=N,
					group.variable=outcome.variable,
					outcome.variable=outcome.variable,
					number.ss.samples.per.iteration=number.ss.samples.per.iteration)
			
			group.names <- sort(unique(as.character(outcome)))
			estimate <- rep(0,length(group.names))
			names(estimate) <- group.names
			for(i in seq(along=group.names)){
				outi <- 1*(outcome==group.names[i])
				
				weightsi <- weights.nomiss
				if(is.null(to.group0.variable)) {
					deg0 <- rep(NA,length(deg))
					deg1 <- rep(NA,length(deg))
					for(k in 1:nrow(rds.data.nomiss)){
						o=outi[rid == id[k]]
						if(length(o)>0){
							deg0[k] = sum(!o) * deg[k] / length(o)
							deg1[k] = sum( o) * deg[k] / length(o)
						}else{
							weightsi[k] <- 0
						}
					}
					weightsi[is.na(deg0)|is.na(deg1)] <- 0
					deg0[is.na(deg0)] <- 0
					deg1[is.na(deg1)] <- 0
				}else{
					deg0 <- as.vector(rds.data.nomiss[,to.group0.variable])
					deg1 <- as.vector(rds.data.nomiss[,to.group1.variable])
				}
				prop.outcome <- HT.estimate(weights=weightsi,outcome=outi)
				mean.degree.outcome0 <- HT.estimate(weights=weightsi,outcome=deg*(outi==0)) / (1-prop.outcome)
				mean.degree.outcome1 <- HT.estimate(weights=weightsi,outcome=deg*(outi==1)) / prop.outcome
				mean.degree <- HT.estimate(weights=weightsi,outcome=deg)
				s01 <- ( outi)*deg0 + (!outi)*deg1
				mean.s01 <- HT.estimate(weights=weightsi,outcome=s01)
				estimate[i] <- 2*prop.outcome*mean.degree.outcome0*(1-prop.outcome)*mean.degree.outcome1 / (mean.s01*mean.degree)
			}
		}
		class(estimate)<-"table" 
	}
	
	result <- homophily.estimate(estimate,outcome.variable,weight.type,recruitment)

	return(result)
}

#' This function computes an estimate of the population homophily and the recruitment homophily based on 
#' a categorical variable.
#' @param rds.data An \code{rds.data.frame} that indicates recruitment patterns by a pair of attributes named ``id'' and ``recruiter.id''.
#' @param outcome.variable A string giving the name of the variable in the \code{rds.data} that contains a categorical or numeric variable to be analyzed.
#' @param weight.type A string giving the type of estimator to use. The options are\cr
#' \code{"Gile's SS"}, \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I/DS"}, \code{"Good-Fellows"} and \code{"Arithemic Mean"}. If \code{NULL} it defaults to \code{"Gile's SS"}.
#' @param uncertainty A string giving the type of uncertainty estimator to use. The options are \code{"Gile's SS"} and \code{"Salganik"}. This is usually determined by \code{weight.type} to be consistent with the estimator's origins (e.g., for \code{"Gile's SS"}, \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I/DS"}, and \code{"Arithemic Mean"}). Hence it's current functionality is limited. If \code{NULL} it defaults to \code{"Gile's SS"}.
#' @param recruitment A logical indicating if the homophily in the recruitment chains should be computed also. The default is FALSE.
#' @param N An estimate of  the number of members of the population being sampled. If \code{NULL} it is read as the \code{population.size.mid} attribute of the \code{rds.data} frame. If that is missing it defaults to 1000.
#' @param to.group0.variable The number in the network of each survey respondent who have group variable value 0. Usually this is not available. The default is to not use this variable.
#' @param to.group1.variable The number in the network of each survey respondent who have group variable value 1. Usually this is not available. The default is to not use this variable.
#' @param number.ss.samples.per.iteration The number of samples to take in estimating the inclusion probabilites in each iteration of the sequential sampling algorithm. If \code{NULL} it is read as the\cr
#' \code{number.ss.samples.per.iteration} attribute of \code{rds.data}. If that is missing it defaults to 5000.
#' @param confidence.level The confidence level for the confidence intervals. The default is 0.95 for 95\%.
#' @return If \code{outcome.variable} is binary then the homophily estimate of
#' 0 verses 1 is returned, otherwise a vector of differential homophily
#' estimates is returned.
#' @section Recruitment Homophily: The recruitment homophily is a homophily
#' measure for the recruitment process. It addresses the question: Do
#' respondents differential recruit people like themselves? That is, the
#' homophily on a variable in the recruitment chains. Take as an example infection
#' status. In this case, it is the ratio of number of recruits that have the
#' same infection status as their recruiter to the number we would expect if there
#' was no homophily on infection status. The difference with the Population
#' Homophily (see below) is that this is in the recruitment chain rather than
#' the population of social ties. For example, of the recruitment homophily on
#' infection status is about 1, we see little effect of recruitment homophily on infection
#' status (as the numbers of homophilous pairs are close to what we would
#' expect by chance).
#' @section Population Homophily: This is an estimate the homophily of a given variable
#' in the underlying networked
#' population. For example, consider HIV status. The
#' population homophily is the homophily in the HIV status of
#' two people who are tied in the underlying population social network (a
#' ``couple''). Specifically, the population homophily is the ratio of the expected
#' number of HIV discordant couples absent homophily to the expected number of HIV
#' discordant couples with the homophily. Hence larger values of population
#' homophily indicate more homophily on HIV status. For example, a value of 1
#' means the couple are random with respect to HIV status. A value of 2 means
#' there are twice as many HIV discordant couples as we would expect if there was
#' no homophily in the population. This measure is meaningful across different
#' levels of differential activity. As we do not see most of the population
#' network, we estimate the population homophily from the RDS data. As an example,
#' suppose the population homophily on HIV is 0.75 so there are 25\% more HIV
#' discordant couples than expected due to chance. So their is actually
#' heterophily on HIV in the population. If the population homophily on sex is
#' 1.1, there are 10\% more same-sex couples than expected due to chance. Hence
#' there is modest homophily on sex.
#' @author Mark S. Handcock with help from Krista J. Gile
#' @references Gile, Krista J., Handcock, Mark S., 2010,
#' \emph{Respondent-driven Sampling: An Assessment of Current Methodology}.
#' Sociological Methodology 40, 285-327.
#' @examples 
#' \dontrun{
#' data(fauxmadrona)
#' names(fauxmadrona)
#' #
#' # True value:
#' #
#' if(require(network)){
#' 	a=as.sociomatrix(fauxmadrona.network)
#' 	deg <- apply(a,1,sum)
#' 	dis <- fauxmadrona.network \%v\% "disease"
#' 	deg1 <- apply(a[dis==1,],1,sum)
#' 	deg0 <- apply(a[dis==0,],1,sum)
#' 	# differential activity
#' 	mean(deg1)/ mean(deg0)
#' 	p=mean(dis)
#' 	N=1000
#' 	# True homophily
#' 	p*(1-p)*mean(deg0)*mean(deg1)*N/(mean(deg)*sum(a[dis==1,dis==0]))
#' }
#' # HT based estimators using the to.group information
#' data(fauxmadrona)
#' homophily.estimates(fauxmadrona,outcome.variable="disease",
#'   to.group0.variable="tonondiseased", to.group1.variable="todiseased",
#'   N=1000)
#' # HT based estimators not using the to.group information
#' homophily.estimates(fauxmadrona,outcome.variable="disease",
#'   N=1000,weight.type="RDS-II")
#' }
#'@export
homophily.estimates <- function(rds.data,outcome.variable,weight.type=NULL,uncertainty=NULL,
		recruitment=FALSE,N=NULL,to.group0.variable=NULL, to.group1.variable=NULL,
		number.ss.samples.per.iteration=NULL,confidence.level=.95)
{
	subset <- substitute(subset)
	if(length(outcome.variable) == 1){
		result <- homophily.estimates.local(rds.data=rds.data,
				outcome.variable=outcome.variable,
				weight.type=weight.type,
				recruitment=recruitment,
				N=N,
				to.group0.variable=to.group0.variable,
				to.group1.variable=to.group1.variable,
				number.ss.samples.per.iteration=number.ss.samples.per.iteration,
				confidence.level=confidence.level)
	}
	else {
		result <- lapply(X=outcome.variable,FUN=function(g){homophily.estimates.local(
							rds.data=rds.data,
							outcome.variable=g,
							weight.type=weight.type,
							recruitment=recruitment,
							N=N,
							to.group0.variable=to.group0.variable,
							to.group1.variable=to.group1.variable,
							number.ss.samples.per.iteration=number.ss.samples.per.iteration,
							confidence.level=confidence.level)})
		names(result) <- outcome.variable
		
	}
	return(result)
}


summary.h.table <- function (object, ...) 
{
	if (!inherits(object, "table")) 
		stop("'object' must inherit from class \"table\"")
	n.cases <- sum(object)
	n.vars <- length(dim(object))
	y <- list(n.vars = n.vars, n.cases = n.cases)
	if (n.vars > 1) {
		m <- vector("list", length = n.vars)
		relFreqs <- object/n.cases
		for (k in 1L:n.vars) m[[k]] <- apply(relFreqs, k, sum)
		expected <- apply(do.call("expand.grid", m), 1L, prod) * 
				n.cases
		statistic <- sum((c(object) - expected)^2/expected, na.rm=TRUE)
		lm <- vapply(m, length, 1L)
		parameter <- prod(lm) - 1L - sum(lm - 1L)
		y <- c(y, list(statistic = statistic, parameter = parameter, 
						approx.ok = all(expected >= 5), p.value = stats::pchisq(statistic, 
								parameter, lower.tail = FALSE), call = attr(object, 
								"call")))
	}
	class(y) <- "summary.table"
	y
}
