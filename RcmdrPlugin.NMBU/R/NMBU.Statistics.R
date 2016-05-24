# $Id: NMBU.Statistics.R 35 2014-01-10 21:17:26Z khliland $

## Statistical functions for RcmdrPlugin.NMBU

##################################
## Methods


#############################
# Linear discriminant scores
DA.coef <- function(){ # Technicality of the R Commander
	doItAndPrint('DA.scores()')
}
DA.scores <- function(object=NULL){
	if(is.null(object)){
		try(eval(parse(text=paste("object <- ", activeModel(), sep=""))))
	}
	if(class(object) != "lda"){
		stop("Only available for linear discriminant analysis (LDA)")
	}
	data.name <- object$call$data
	var.names <- colnames(object$means) # Find and choose correct data
	X <- as.matrix(eval(parse(text = paste(data.name, "[,c('", paste(var.names, collapse="','"), "')]", sep=""))))
	n <- dim(X)

	# Initialisations
	y       <- eval(parse(text=paste(ActiveDataSet(), "$",formula(object)[[2]])))
	groups  <- levels(y)
	g       <- length(groups)
	n.g     <- numeric(g)
	index.g <- list()
	Y.dummy <- matrix(0,n[1],g)
	for(i in 1:g){
		Y.dummy[,i]  <- y==groups[i]
		n.g[i]       <- sum(Y.dummy[,i])
		index.g[[i]] <- which(Y.dummy[,i]==1)
	}
	colnames(Y.dummy) <- groups
	names(n.g) <- groups

	the.means <- matrix(0, n[2], g)
	the.covs  <- list()
	for(i in 1:g){ # Compute group-wise means and covariances
		the.means[,i] <- apply(X[index.g[[i]],,drop=FALSE],2,mean)
		the.covs[[i]] <- var(X[index.g[[i]],])
	}
	S <- matrix(0, n[2],n[2])
	for(i in 1:g){ # Common covariance matrix
		S <- S + (n.g[i]-1)/(n[1]-g)*the.covs[[i]]
	}
	
# Common classification for two classes
	if(g == 2 && class(object) == "lda"){
		beta0 <- -0.5*t(the.means[,1]-the.means[,2]) %*% solve(S) %*% (the.means[,1]+the.means[,2]) - log(object$prior[2]/object$prior[1])
		beta  <- t(the.means[,1]-the.means[,2]) %*% solve(S)
		cat('Classify as ', groups[1], ' if \n', beta0, sep='')
		for(j in 1:n[2]){
			cat(' + ', beta[j], '*', var.names[j], sep='')
		}
		cat(' > 0\n\n', sep='')
	}
	
# Group-wise scores
	cat("Group-wise linear discriminant scores")
	if(class(object) == "lda"){
		for(i in 1:g){
			beta0 <- log(object$prior[i]) - 0.5*t(the.means[,i])%*%solve(S)%*%the.means[,i]
			beta  <- the.means[,i] %*% solve(S)
			cat("\nd.", groups[i], " = ", beta0, sep="")
			for(j in 1:n[2]){
				cat(" + ", beta[j], "*", var.names[j], sep="")
			}
		}
	}
	cat("\n")
}


#################################
# Hierarchical clustering merges
hclust.list <- function(){ # Stupid technicality of the R Commander
	doItAndPrint('hclust.merge()')
}
hclust.merge <- function(object){
	if(missing(object)){
		solutionNumber<-length(listHclustSolutions())
		try(eval(parse(text=paste("object <- HClust.", solutionNumber, sep=""))))
	}
	cat(paste("Linkage: ", object$method, "\n", sep=""))
	print(data.frame(Merge=object$merge, Height=object$height))
}


################################
# Prediction for lm CI/PI
predict_CI_PI <- function(model, data, level, xXXx=FALSE){
	.activeDataSet <- ActiveDataSet()
	data.names <- colnames(data)
	for(i in 1:length(data.names)){
		if(justDoIt(paste("is.factor(",.activeDataSet,"$", data.names[i], ")", sep=""))&& is.numeric(data[,i])){
			data[,i] <- as.factor(data[,i])
		}
	}
	CI <- predict(model, data, interval='conf', level=level)
	PI <- predict(model, data, interval='predict', level=level)
	CIPI <- cbind(CI,PI[,2:3, drop=FALSE])
	colnames(CIPI) <- c('fit','lwr.CI','upr.CI','lwr.PI','upr.PI')
	if(xXXx){
		the.model <- justDoIt(ActiveModel())
		X <- model.matrix(the.model)
		x <- model.matrix(delete.response(terms(the.model)),model.frame(delete.response(terms(the.model)),data,xlev=the.model$xlevels))
		n <- dim(x)[1]
		result <- numeric(n)
		for(i in 1:n){
			result[i] <- x[i,] %*% solve(t(X) %*% X) %*% x[i,]
		}
		CIPI <- cbind(CIPI,h_00=result)

	}
	CIPI
}

################################
# Prediction for glm link/response
predict_link_response <- function(model, data){
	.activeDataSet <- ActiveDataSet()
	data.names <- colnames(data)
	for(i in 1:length(data.names)){
		if(justDoIt(paste("is.factor(",.activeDataSet,"$", data.names[i], ")", sep=""))&& is.numeric(data[,i])){
			data[,i] <- as.factor(data[,i])
		}
	}
	link     <- predict(model, data, type='link', se.fit=TRUE)
	response <- predict(model, data, type='response', se.fit=TRUE)
	LS <- cbind(link$fit, link$se.fit, response$fit, response$se.fit)
	colnames(LS) <- c('link','se.link','response','se.response')
	LS
}


################################
# Extended numerical summary
sumSq <- function(x,na.rm = TRUE){
	if(is.vector(x))
		sum(x^2,na.rm = TRUE)
	else
		apply(x^2,2,sum, na.rm=TRUE)
}
sdErr <- function(x,na.rm = TRUE){
	if(is.vector(x))
		sd(x,na.rm = TRUE)/sqrt(sum(!is.na(x)))
	else
		apply(x,2,sd,na.rm = TRUE)/sqrt(apply(!is.na(x), 2, sum, na.rm=TRUE))
}
numSummaryNMBU <- function(data, 
		statistics=c("mean", "median", "sum", "sumSq", "sd", "sdErr", "var", "quantiles", "cv", "skewness", "kurtosis"),
		type=c("1", "2", "3"),
		quantiles=c(0, .25, .5, .75, 1), groups){
	sd <- function(x, type, ...){
		apply(as.matrix(x), 2, stats::sd, na.rm=TRUE)
	}
	cv <- function(x, ...){
		mean <- mean(x, na.rm=TRUE)
		sd <- sd(x)
		if (any(x <= 0, na.rm=TRUE)) warning("not all values are positive")
		cv <- sd/mean
		cv[mean <= 0] <- NA
		cv
	}
	skewness <- function(x, type, ...){
		if (is.vector(x)) return(e1071::skewness(x, type=type, na.rm=TRUE))
		apply(x, 2, skewness, type=type)
	}
	kurtosis <- function(x, type, ...){
		if (is.vector(x)) return(e1071::kurtosis(x, type=type, na.rm=TRUE))
		apply(x, 2, kurtosis, type=type)
	}
	data <- as.data.frame(data)
	if (!missing(groups)) groups <- as.factor(groups)
	variables <- names(data)
	if (missing(statistics)) statistics <- c("mean", "sd", "quantiles")
	statistics <- match.arg(statistics, c("mean", "median", "sum", "sumSq", "sd", "sdErr", "var", "quantiles", "cv", "skewness", "kurtosis"),
			several.ok=TRUE)
	type <- match.arg(type)
	type <- as.numeric(type)
	ngroups <- if(missing(groups)) 1 else length(grps <- levels(groups))
	quantiles <- if ("quantiles" %in% statistics) quantiles else NULL
	quants <- if (length(quantiles) > 1) paste(100*quantiles, "%", sep="")
			else NULL
	nquants <- length(quants)
	stats <- c(c("mean", "median", "sum", "sumSq", "sd", "sdErr", "var", "cv", "skewness", "kurtosis")[c("mean", "median", "sum", "sumSq", "sd", "sdErr", "var", "cv", "skewness", "kurtosis") %in% statistics], quants)
	nstats <- length(stats)
	nvars <- length(variables)
	result <- list()
	if ((ngroups == 1) && (nvars == 1) && (length(statistics) == 1)){
		if (statistics == "quantiles"){
			table <- quantile(data[,variables], probs=quantiles, na.rm=TRUE,type=type)
		} else {
			if (statistics == "kurtosis" || statistics == "skewness"){
				table <- do.call(statistics, list(x=data[,variables], na.rm=TRUE, type=type))
			} else {
				table <- do.call(statistics, list(x=data[,variables], na.rm=TRUE))
			}
			names(table) <- statistics
		}
		NAs <- sum(is.na(data[,variables]))
		n <- nrow(data) - NAs
		result$type <- 1
	}
	else if ((ngroups > 1)  && (nvars == 1) && (length(statistics) == 1)){
		if (statistics == "quantiles"){
			table <- matrix(unlist(tapply(data[, variables], groups,
									quantile, probs=quantiles, na.rm=TRUE,type=type)), ngroups, nquants,
					byrow=TRUE)
			rownames(table) <- grps
			colnames(table) <- quants
		} else {
			if( statistics == "skewness" || statistics == "kurtosis" ){
				table <- tapply(data[,variables], groups, statistics,
					na.rm=TRUE, type=type)
			}
				table <- tapply(data[,variables], groups, statistics,
					na.rm=TRUE)
		}
		NAs <- tapply(data[, variables], groups, function(x)
					sum(is.na(x)))
		n <- table(groups) - NAs
		result$type <- 2
	}
	else if ((ngroups == 1) ){
		X <- as.matrix(data[, variables])
		table <- matrix(0, nvars, nstats)
		rownames(table) <- if (length(variables) > 1) variables else ""
		colnames(table) <- stats
		if ("mean" %in% stats) table[,"mean"] <- apply(X, 2, mean, na.rm=TRUE)
		if ("median" %in% stats) table[,"median"] <- apply(X, 2, median, na.rm=TRUE)
		if ("sum" %in% stats) table[,"sum"] <- apply(X, 2, sum, na.rm=TRUE)
		if ("sumSq" %in% stats) table[,"sumSq"] <- sumSq(X, na.rm=TRUE)
		if ("sd" %in% stats) table[,"sd"] <- apply(X, 2, sd, na.rm=TRUE)
		if ("sdErr" %in% stats) table[,"sdErr"] <- sdErr(X, na.rm=TRUE)
		if ("var" %in% stats) table[,"var"] <- apply(X, 2, var, na.rm=TRUE)
		if ("cv" %in% stats) table[,"cv"] <- cv(X)
		if ("skewness" %in% statistics) table[, "skewness"] <- skewness(X, type=type)
		if ("kurtosis" %in% statistics) table[, "kurtosis"] <- kurtosis(X, type=type)
		if ("quantiles" %in% statistics){
			table[,quants] <- t(apply(data[, variables, drop=FALSE], 2, quantile,
							probs=quantiles, na.rm=TRUE,type=type))
		}
		NAs <- colSums(is.na(data[,variables, drop=FALSE]))
		n <- nrow(data) - NAs
		result$type <- 3
	}
	else {
		table <- array(0, c(ngroups, nstats, nvars),
				dimnames=list(Group=grps, Statistic=stats, Variable=variables))
		NAs <- matrix(0, nvars, ngroups)
		rownames(NAs) <- variables
		colnames(NAs) <- grps
		for (variable in variables){
			if ("mean" %in% stats)
				table[, "mean", variable] <- tapply(data[, variable],
						groups, mean, na.rm=TRUE)
			if ("median" %in% stats)
				table[, "median", variable] <- tapply(data[, variable],
						groups, median, na.rm=TRUE)
			if ("sum" %in% stats)
				table[, "sum", variable] <- tapply(data[, variable],
						groups, sum, na.rm=TRUE)
			if ("sumSq" %in% stats)
				table[, "sumSq", variable] <- tapply(data[, variable],
						groups, sumSq, na.rm=TRUE)
			if ("sd" %in% stats)
				table[, "sd", variable] <- tapply(data[, variable],
						groups, sd, na.rm=TRUE)
			if ("sdErr" %in% stats)
				table[, "sdErr", variable] <- tapply(data[, variable],
						groups, sdErr, na.rm=TRUE)
			if ("var" %in% stats)
				table[, "var", variable] <- tapply(data[, variable],
						groups, var, na.rm=TRUE)
			if ("cv" %in% stats)
				table[, "cv", variable] <- tapply(data[, variable],
						groups, cv)
			if ("skewness" %in% stats)
				table[, "skewness", variable] <- tapply(data[, variable],
						groups, skewness, type=type)
			if ("kurtosis" %in% stats)
				table[, "kurtosis", variable] <- tapply(data[, variable],
						groups, kurtosis, type=type)
			if ("quantiles" %in% statistics) {
				res <- matrix(unlist(tapply(data[, variable], groups,
										quantile, probs=quantiles, na.rm=TRUE, type=type)), ngroups, nquants,
						byrow=TRUE)
				table[, quants, variable] <- res
			}
			NAs[variable,] <- tapply(data[, variable], groups, function(x)
						sum(is.na(x)))
		}
		if (nstats == 1) table <- table[,1,]
		if (nvars == 1) table <- table[,,1]
		n <- table(groups)
		n <- matrix(n, nrow=nrow(NAs), ncol=ncol(NAs), byrow=TRUE)
		n <- n - NAs
		result$type <- 4
	}
	result$table <- table
	result$statistics <- statistics
	result$n <- n
	if (any(NAs > 0)) result$NAs <- NAs
	class(result) <- c("numSummary")
	result
}
xtable.numSummary <- function(x, caption = NULL, label = NULL, align = NULL, digits = NULL, 
                              display = NULL, ...){
  xtable(x$table, caption, label, align, digits, display, ...)
}


################################
# Goodness of fit (multinom and polr)
deviance_tests <- function(){
	tmp.mod <- eval(parse(text=ActiveModel()))
	if(formula(tmp.mod)[[3]] != 1){
		if(any(class(tmp.mod)=="polr")){
			if(any(names(tmp.mod$call)=="weights")){
				saturated.m <- eval(parse(text = paste("multinom(", formula(tmp.mod)[[2]], "~", formula(tmp.mod)[[3]], ", weights=", tmp.mod$call['weights'][[1]], ", data=", ActiveDataSet(), ",Hess=TRUE)", sep="")))
			} else {
				saturated.m <- multinom(formula(tmp.mod),data=eval(parse(text=ActiveDataSet())),Hess=TRUE)
			}
			saturated.m$df.residual <- tmp.mod$n-saturated.m$edf
			class(saturated.m) <- "polr"
		} else {
			saturated.m <- update(tmp.mod,paste("~ factor(",paste(fparse(formula(tmp.mod)),sep="",collapse=")*factor("),")",sep=""))
		}
		intercept.m <- update(tmp.mod,~1)
		if(is.logical(all.equal(saturated.m,tmp.mod))){
			A <- anova(intercept.m, tmp.mod, test="Chisq")
			if(colnames(A)[1]=="Model")
				A[1] <- c("No effects",ActiveModel())
		} else {
			A <- anova(intercept.m, tmp.mod, saturated.m, test="Chisq")
			if(colnames(A)[1]=="Model")
				A[1] <- c("No effects",ActiveModel(),"Saturated")
		}
	} else {
		A <- NULL
	}
	A
}


################################
# Extended summary from multinom
summaryMultinom <- function(object){
	M <- summary(object, Wald=TRUE)
	cat('Call:\n')
	print(M$call)
	if(terms(object)[[3]]!=1){
		A <- Anova(object)}
	N <- dim(M$Wald.ratios)
	n <- prod(N)
	lH <- numeric(n)
	zeros <- rep(0,n)
	if(!is.null(N) && n>1){
		for(i in 1:n){
			z.tmp <- zeros
			z.tmp[i] <- 1
			lH[i] <- linearHypothesis(object, z.tmp, rhs=0)[2,3]
		}
		p.values <- matrix(lH,N[1],N[2],byrow=TRUE)
		result <- cbind(as.vector(t(M$coefficients)),as.vector(t(M$standard.errors)),as.vector(t(M$Wald.ratios)),as.vector(t(p.values)))
		colnames(result) <- c('Coef','SE Coef','Z','P')
	} else {
		result <- cbind(as.vector(t(M$coefficients)),as.vector(t(M$standard.errors)),as.vector(t(M$Wald.ratios)))
		colnames(result) <- c('Coef','SE Coef','Z')
	}
	rn <- rownames(M$coefficients)
	cn <- colnames(M$coefficients)
	j <- 0
	if(!is.null(N) && n>1){
		for(i in 1:N[1]){
			R <- result[1:N[2]+j,,drop=FALSE]
			dimnames(R) <- eval(parse(text=paste("list('", rn[i], "'=c('", paste(cn,collapse="','",sep=""), "'),", "c('Coef','SE Coef','Z','P'))", sep="")))
			j <- j+N[2]
			print(R)
		}
	}
	cat("\n")
	if(terms(object)[[3]]!=1){
		print(A)}
	cat(paste("\nResidual Deviance: ", round(M$deviance,4), "\n", sep=""))
	cat(paste("AIC: ", round(M$AIC,4), "\n", sep=""))
	print(logLik(object))
}

################################
# Extended summary from multinom
extend.colnames <- function(object, the.name){
	if(is.matrix(object)){
		colnames(object) <- paste(the.name,colnames(object),sep='.')
	}
	object
}


################################
# Extended summary polr
summaryOrdinal <- function(object, digits = max(3, .Options$digits - 3),
                         correlation = FALSE, ...)
{
	if(terms(object)[[3]]!=1){
		A <- Anova(object)}
	lL <- logLik(object)
    cc <- c(coef(object), object$zeta)
    pc <- length(coef(object))
    q <- length(object$zeta)
    coef <- matrix(0, pc+q, 4L, dimnames=list(names(cc),
                               c("  Coef", " SE Coef", "  Z", "  P")))
    coef[, 1L] <- cc
    vc <- vcov(object)
    coef[, 2L] <- sd <- sqrt(diag(vc))
    coef[, 3L] <- coef[, 1L]/coef[, 2L]
	coef[, 4L] <- 1-pchisq(coef[, 3L]^2,1)
    object$coefficients <- coef
    object$pc <- pc
    object$digits <- digits
    if(correlation)
        object$correlation <- (vc/sd)/rep(sd, rep(pc+q, pc+q))
    class(object) <- "summary.polr"
    print(object)
	print(lL)
	cat("\n")
	if(terms(object)[[3]]!=1){
		print(A)}
}



