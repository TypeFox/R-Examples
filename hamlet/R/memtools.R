#
# Tools for the analysis of mixed-effects models (MEM) from e.g. lme4, nlme
#

# Plot residuals of a mixed-effects model along with linear trend fit
mem.plotresid <- function(fit, linear=T, type="XbZu", main, xlab, ylab){
	res <- mem.getcomp(fit)
	plot(res[,type], res[,"e"], pch=16, xlab="", ylab="")
	
	if(linear){
		residlm <- lm(e ~ y, data.frame(y = res[,type], e = res[,"e"]))
		abline(residlm)	
		legend("bottomright", lwd=1, col="black", legend=c("Linear trend"), bty="n")	
	}
	if(!missing(main)) title(main=main)
	if(!missing(xlab)){ title(xlab=xlab)}
	else{ title(xlab=type)}
	if(!missing(ylab)){ title(ylab=ylab)}
	else{ title(ylab="e")}
}

# Plot histogram distributions of random effects
mem.plotran <- function(fit, breaks=100){
	count <- 0
	if(class(fit) %in% c("lmerMod", "mer", "merModLmerTest", "nlmerMod")){
		lapply(lme4::ranef(fit), FUN=function(x) count <<- count + ncol(x))
		par(mfrow=c(ceiling(sqrt(count)),ceiling(sqrt(count))))
		lapply(lme4::ranef(fit), FUN=function(x) apply(x, MARGIN=2, FUN=function(y) hist(y, breaks=breaks)))
	}else if(class(fit) %in% c("nlme", "lme", "merModLmerTest")){
		count <- count + ncol(nlme::ranef(fit))
		par(mfrow=c(ceiling(sqrt(count)),ceiling(sqrt(count))))
		apply(nlme::ranef(fit), MARGIN=2, FUN=function(y) hist(y, breaks=breaks))
	}else{
		stop(paste("Invalid class of fit:", class(fit)))
	}
	
}

# Get per-observation components Xb and Zu from a mixed-effects model
mem.getcomp <- function(fit){
	if(class(fit) %in% c("lmerMod", "glmerMod", "mer", "merModLmerTest")){
		# Xb 
		fix <- lme4::getME(fit,'X') %*% lme4::fixef(fit)
		# Zu
		ran <- t(as.matrix(lme4::getME(fit,'Zt'))) %*% unlist(lme4::ranef(fit))
		# Xb + Zu
		fixran <- lme4::getME(fit,'y') - resid(fit)
		# y
		y = lme4::getME(fit, 'y')
	}else if(class(fit)=="nlmerMod"){ # Non-linear lme4-fit
		# Xb 
		fix <- lme4::getME(fit,'X') %*% lme4::fixef(fit)
		# Zu
		ran <- t(as.matrix(lme4::getME(fit,'Zt'))) %*% unlist(lme4::ranef(fit))
		# Xb + Zu
		fixran <- lme4::getME(fit,'y') - resid(fit)
		# y
		y = lme4::getME(fit, 'y')
	}else if(class(fit)=="lme"){
		# Xb
		fix <- fit$fitted[,1]
		# Xb + Zu
		fixran <- fit$fitted[,2]
		# Zu = Xb + Zu - Xb
		ran <- fixran - fix
		# y
		y <- fit$data[,as.character(fit$terms[[2]])]
	}else if(class(fit)=="lm"){
		fix = fixran = fitted(fit)
		ran = 0
		y <- fit$model[,as.character(fit$terms[[2]])]
	}else{
		stop(paste("Invalid class of fit:", class(fit)))
	}
	result <- cbind(fix = fix, ran = ran, fixran = fixran, e = resid(fit), y = y)
	colnames(result) <- c("Xb", "Zu", "XbZu", "e", "y")
	result
}

# Bootstrap sampled data from a fitted lme4-object and compute power curves for the corresponding fixed-effects terms by re-fitting the same model
mem.powersimu <- function(
	# Fitted lme4 object
	fit,
	
	# N-values to be tested using bootstrap
	# Corresponds to number of unique individuals (usually extracted from random effects formulation)
	N=4:20,
	# Number of bootstraps per each N
	boot = 100,
	
	# Factor to use as a grouping variable to sample from; if NULL the first identified factor from the lme4-object will be used
	# This should correspond to a feasible column name in the original data frame provided to fit the data, i.e. unique label column
	level = NULL,
	# If strata should be adjusted; for non-paired models, this usually should be the treatment groups if they are desired to be of equal size
	# Should correspond to a variable available in the model data frame; each unique level will be sampled in equal amounts
	strata = NULL,
	
	# If a model cannot be fitted, what is the default inference for coefficient significance; 
	# FALSE means in non-informative cases, a conservative estimate favoring non-significance is given, 
	# TRUE means an overestimated significance is given
	default = FALSE,
	# Set a random seed for bootstrapping
	seed = NULL,
	
	# Should a power-curve be plotted based on the power matrix
	plot = TRUE,
	# Should a loess-smoothed power-curve be plotted alongside the actual observed sampled powers
	plot.loess = FALSE,
	# If a legend is plotted, a text field indicating where it should be placed in the figure
	legendpos = "bottomright",
	
	# Should bootstrapped data be returned;
	# Notice that this will be _highly_ inefficient if large amounts of bootstrapping is performed
	# It is however useful for inspecting and making sure the bootstrapping schema is working correctly in the context of the fitted lme4-model
	return.data = FALSE,
	# Level of verbosity, 0: none; 1: normal; 2: debugging
	verb = 1,
	# Additional parameters
	...
){
	if(!class(fit) %in% c("lme4", "lmerMod", "merModLmerTest")){
		stop(paste("Invalid class of fit (should be a lme4-object), given class:", class(fit)))
	}
	# Set random seed if desired
	if(!is.null(seed)) set.seed(seed)
	
	# If used desires bootstrapped returned data, it will be returned using lists of lists with bootstrapped data frames
	if(return.data) dat <- list()
	
	# Extract model formula from the lme4-object	
	form = attr(fit@frame, "formula")
	# Extract data frame
	frame = fit@frame
	# Levels of factors; used to determine the sampling unit
	factors = fit@flist
	# Identify the sampling unit, usually the identification code for individuals
	if(is.null(level)){
		level = names(factors)[1]
		if(!length(names(factors))==1) warning(paste("Multiple identification codes present, but using only", level))
	}
	# Identifying sampling strata, for example intervention groups that should be balanced in respect to the amount of individuals
	if(!is.null(strata)){
		if(!strata %in% colnames(frame)){
			stop(paste("Stratifying coefficient",strata,"is not present in the model frame"))
		}else{
			tab = table(frame[,c(strata, level)])
			ids = lapply(1:nrow(tab), FUN=function(z) colnames(tab)[tab[z,]>0])
		}
	}else{
		# Unique IDs to sample from
		#ids = levels(factors[[level]])
		ids = list(levels(factors[[level]]))
	}
	# Sequence of total counts of individuals;
	# for example if doing balanced sampling of 2, 5 and 10 individuals from 2 substrata,
	# the total N coulds will be 4, 10 and 20 respectively
	Ntot <- N*length(ids)
	if(verb>=2) print(ids)

	# If user does not just wish to obtain returned bootstrapped datasets, re-fitted lme4 models are generated and fixed effects inference is returned	
	if(!return.data){
		# Start the sampling process
		# Loop over N values	
		npmatrix <- do.call("rbind", lapply(N, FUN=function(n){
			# Loop over different bootstrapped samplesets
			p = do.call("cbind", lapply(1:boot, FUN=function(x){
				# Bootstrapped identification codes
				#s = sample(ids, size=n, replace=T)
				# Sample so that strata are evenly balanced
				s = unlist(lapply(1:length(ids), FUN=function(z) sample(ids[[z]], size=n, replace=T)))
				# Make new unique identification codes; repeated IDs will be separated by suffix .1, .2, .3, ...
				suniq <- make.unique(as.character(s))
				# Generate new bootstrapped sampleset
				sampleset <- do.call("rbind", lapply(1:length(s), FUN=function(z){
					# Pick sampled individuals
					tmp <- frame[which(s[z]==factors[[level]]),]
					# Original individual identification codes have been make.unique'd due to sampling with replacement
					tmp[,level] <- suniq[z]
					# Return sampled observations and bind to a new data frame
					tmp
				}))

				# Try re-fit the lme4 using same the model formula with a new bootstrapped dataset
				if(verb>=2) print(sampleset)
				newfit <- try(lme4::lmer(form, data = sampleset), silent=T)
				# If non-error, presumably a feasible fit
				if(!class(newfit) %in% c("try-error")){
					# Extract p-values using the Satterthwaite approximation in lmerTest-package
					summarized <- lmerTest::summary(newfit)$coefficients
					# If p-values could be computed in lmerTest
					if("Pr(>|t|)" %in% colnames(summarized)){
						summarized[,"Pr(>|t|)"]<0.05
					# Else if, test statistic approximation using significant: |t|>2
					}else if("t value" %in% colnames(summarized)){
						abs(summarized[,"t value"])>2
					# Else a fitting error, utilize default inference
					}else{
						rep(default, times=length(lme4::fixef(fit)))
					}
				}
				# Else a fitting error, utilize default inference
				else{
					rep(default, times=length(lme4::fixef(fit)))
				}
			}))
			apply(p, MARGIN=1, FUN=function(z){
				sum(z)/length(z)
			})
		}))
		# Name the rows and columns of the power simulation matrix
		rownames(npmatrix) <- paste("GroupN_", N, "_TotalN_",Ntot, sep="")
		colnames(npmatrix) <- names(lme4::fixef(fit))
	}
	# Else user only wishes to generate bootstrapped datasets without the re-fitting of the lme4 model
	else{
		# Start the sampling process
		# Loop over N values	
		# In returning data, the outer loop will be just a list over N values
		npmatrix <- lapply(N, FUN=function(n){
			# Loop over different bootstrapped samplesets
			# In returning data, the inner loop will be just the different bootstrapped sets for a specific N value
			tmp <- lapply(1:boot, FUN=function(x){
				# Bootstrapped identification codes
				#s = sample(ids, size=n, replace=T)
				# Sample so that strata are evenly balanced
				s = unlist(lapply(1:length(ids), FUN=function(z) sample(ids[[z]], size=n, replace=T)))
				# Make new unique identification codes; repeated IDs will be separated by suffix .1, .2, .3, ...
				suniq <- make.unique(as.character(s))
				# Generate new bootstrapped sampleset
				sampleset <- do.call("rbind", lapply(1:length(s), FUN=function(z){
					# Pick sampled individuals
					tmp <- frame[which(s[z]==factors[[level]]),]
					# Original individual identification codes have been make.unique'd due to sampling with replacement
					tmp[,level] <- suniq[z]
					# Return sampled observations and bind to a new data frame
					tmp
				}))
				sampleset
			})
			names(tmp) <- paste("boot_", 1:boot, sep="")
			tmp
		})
		names(npmatrix) <- paste("N_", N, sep="")
		npmatrix
	}
	
	# If default power curve plotting should be done
	if(plot & !return.data){
		plot.new()
		#plot.window(xlim=range(N), ylim=range(npmatrix))
		plot.window(xlim=range(Ntot), ylim=c(0,1))
		box(); axis(1); axis(2); title(xlab="Total N", ylab="Power")
		# Plot power points and the corresponding curve
		for(i in 1:ncol(npmatrix)){
			x = Ntot
			y = npmatrix[,i]
			points(x=x, y=y, type="l", col=i)
			points(x=x, y=y, pch=16, col=i)
			# If loess smoothed power curve should be plotted along the actual observed sampled powers
			if(plot.loess){
				xseq = seq(from=Ntot[1], to=Ntot[length(Ntot)], length.out=100)
				l = loess(y~x, data = data.frame(x,y))
				ypred = predict(l, newdata = data.frame(x=xseq))
				points(x=xseq, y=ypred, type="l", lty=2, col=i)
			}
		}
		# Plot legend
		legend(legendpos, bty="n", lwd=1, col=1:ncol(npmatrix), legend=colnames(npmatrix))
		# Loess legend
		if(plot.loess) legend("topright", bty="n", lwd=1, lty=2, legend="Loess smoothed curve")
	}

	# Return the obtained inference
	# For !return.data (default choice), it is the statistical power per each fixed effects coefficient (columns) per N values (rows)
	# for return.data==TRUE, it is the bootstrapped datasets in a list of lists, where outer loop is the N values, and inner loop is the bootstraps
	npmatrix
}

