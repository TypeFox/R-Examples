## AUXILIARY FUNCTIONS

# setNumericPredictors() defines the values of numeric predictors (covariates and offset)
setNumericPredictors <- function(model, numeric.predictors, covariates=NULL, offset=NULL){
	mf <- getModelFrame(model)
	if (length(numeric.predictors) > 0){
		# If no covariates were specified by the user, use the default (average)
		if (is.null(covariates)){
			covariates <- colMeans(mf[numeric.predictors])
		}else{
			# If covariates have names, complete with default (ignoring covariates that are not in the model)
			covariate.names <- names(covariates)
			if (!is.null(covariate.names)){
				valid.covariates <- covariate.names %in% numeric.predictors
				if (!all(valid.covariates)) warning("Some covariates are not in the model, and will be ignored.")
				user.covariates <- covariates[valid.covariates]
				covariates <- colMeans(mf[numeric.predictors])
				covariates[names(user.covariates)] <- user.covariates
			}else{
				# If they are unnamed, assign them in the default order of the model, after adjusting the length of the vector
				covariates = rep(covariates,length=length(numeric.predictors))
				names(covariates) <- numeric.predictors
			}
		}
	} else covariates <- NULL
	# Offset (if any)
	offsetvar <- getOffset(model)
	if (length(offsetvar) > 0){
		# If no offset value is specified by the user, use the default (average)
		offsetvalue <- if(is.null(offset)) mean(offsetvar) else offset
	}else{
		offsetvalue <- 0
		if (!is.null(offset)) warning("The model does not contain an offset. Offset argument will be ignored.")
	}
	list(covariates=covariates,offsetvalue=offsetvalue)
}

# addPredictorsToFrame() is used to include covariates in a data frame
# of factor combinations. Add offset value, if any
addPredictorsToFrame <- function(factor.frame, numeric.predictors, term){
	if (length(numeric.predictors) > 0){
		# Repeat the numeric values along the rows of the between-subjects data frame
		new.columns <- t(matrix(term$covariates,dimnames=list(numeric.predictors)))
		if (nrow(factor.frame)>0){
			as.data.frame(cbind(factor.frame,new.columns))
		}else{
			as.data.frame(new.columns)
		}
	}else{factor.frame}
}

# defineLHT() makes the first approximation to the Linear Hypothesis Matrix
# If the model has predictors, define L according to the model formula
defineLHT <- function(model,term,factor.frame){
	tm <- terms(model)
	if (nrow(factor.frame)){
		# rhf <- formula(terms(model))[c(1,3)] # 0.1-3 definition (included offset)
		# Define right-hand side of the model formula,
		# with masked names of predictors to avoid "AsIs", etc.
		predictors <- names(factor.frame)
		masked.predictors <- make.names(predictors, unique=TRUE)
		masked.terms <- attr(tm, "term.labels")
		# Masking
		for (p in which(predictors != masked.predictors)){
			masked.terms <- gsub(predictors[p], masked.predictors[p], masked.terms, fixed=TRUE)
		}
		rhf <- paste(c(attr(tm,"intercept"), masked.terms), collapse= "+")
		rhf <- as.formula(paste("~", rhf))
		names(factor.frame) <- masked.predictors
		L <- model.matrix(rhf, data=factor.frame)
		colnames(L) <- rownames(as.matrix(getCoef(model)))
	# Otherwise, L is defined to just get the intercept
	}else{
		if (!attr(tm,"intercept")) stop("Null model (no predictors or intercept).")
		L <- matrix(1,dimnames=list(NULL,"(Intercept)"))
	}
	# If the term to analyse is not the intercept,
	# set to zero the rows of the model matrix unrelated to specified covariates
	if (term$num.vars[1] != "(Intercept)"){
		# This matrix tells what terms are affected by each variable
		# We choose the rows for the variables in the selected term
		terms.matrix <- attr(tm,"factors")[term$num.vars,,drop=F]
		# Only the columns where all variables are present are of interest to us
		affected.terms <- which(apply(terms.matrix,2,"all"))
		# The attribute "assign" of the model matrix tells what coefficents are related to each term,
		# and we use affected.terms to select the relevant coefficients
		relevant.coefficients <- attr(getModelMatrix(model),"assign") %in% affected.terms
		# Now set irrelevant rows to zero
		L[,!relevant.coefficients] <- 0
	}
	return(L)
}

# checkFactors() is used to check consistency of level combinations in the
# factors defined in levels, and convert them to numeric matrices
checkFactors <- function(frame,factor.names,levels){
	for (fname in factor.names){
		f <- frame[[fname]]
		combination <- levels[[fname]]
		# Literal expressions must represent one level or a pair of levels 
		if (is.character(combination)){
			combined.levels <- combination
			if (length(combination)==1){
				# If there is one level, evaluate the values at that level
				combination <- 1
				names(combination) <- combined.levels
			}else if (length(combined.levels)==2){
				# If there is a pair of levels, evaluate the contrast
				combination <- c(1,-1)
				names(combination) <- combined.levels
			}else{
				stop("Incorrect number of literal levels for factor \"",fname,"\" (must be 1 or 2).")
			}
		}
		# Check the numeric coefficients of the factor levels
		if (is.numeric(combination)){
			# Make combination into a matrix (in case it was a vector)
			combination <- as.matrix(combination)
			# Unnamed coefficient matrices must have the same rows as levels in the factor
			if (is.null(rownames(combination))){
				if (nrow(combination) != nlevels(f)) stop("Mismatch in the number of unnamed factor levels for \"",fname,"\".")
				rownames(combination) <- levels(f)
			}else{
				# Named coefficients are assigned to a matrix of factor levels, filled in with zeros
				factor.levels <- match(rownames(combination),levels(f))
				combination.copy <- combination
				combination <- matrix(0,nrow=nlevels(f),ncol=ncol(combination.copy))
				if (any(is.na(factor.levels))) stop("Mismatch in the names of factor levels for \"",fname,"\".")
				combination[factor.levels,] <- combination.copy
				rownames(combination) <- levels(f)
			}
		}else{
			stop("Invalid value assigned to factor \"",fname,"\".")
		}
		levels[[fname]] <- combination
		if (is.null(colnames(levels[[fname]]))) colnames(levels[[fname]]) <- paste(fname,as.character(1:ncol(levels[[fname]])),sep="")
	}
	return(levels)
}

# makeT() creates a transformation matrix, used to define the Linear Hypothesis
# matrix (for between-subjects factors) or the response transformation matrix
# (for within-subjects factors), depending on the combinations of factor levels.
# The transformation matrix is created progressively, by "translating" the combinations
# of each factor into matrices that are sequentially multiplied.
makeT <- function(frame,factor.names,levels){
	# First matrix, defined from the identic rows of the between/within-subjects
	# model data frame, when unspecified factors are removed.
	# A dummy column with ones is added to the model data frame, to avoid
	# problems when all columns be eventually removed.
	# (The name of the dummy column is coerced to be different from any other one.)
	dummyname <- paste("z",max(factor.names),sep=".")
	frame[[dummyname]] <- 1
	# All factors are interacting, the transformation matrix (Tm) in the first step is diagonal.
	if (length(factor.names)==ncol(frame)-1){
		m <- nrow(frame)
		Tm <- diag(m)
	}else{
		# Remove columns of unspecified factors
		frame <- frame[,c(factor.names,dummyname)]
		# Vector with a different value for each combination of factors
		frame_1 <- apply(frame,1,paste,collapse="")
		n <- nrow(frame)
		# Subset of unique elements
		frame <- unique(frame)
		frame_1.m <- unique(frame_1)
		m <- nrow(frame)
		# Matrix with coefficients for averaging identical combinations of factors
		# Rows in Mo are the original (repeated) combinations
		To <- matrix(rep(frame_1,each=m),nrow=m)
		# Columns in Mb are the unique combinations
		Tu <- matrix(rep(frame_1.m,n),nrow=m)
		Tm <- (To==Tu) * m/n
	}
	# Number of contrasts to calculate for the current factor (initialized as 1)
	nc <- 1
	# Labels for the final transformation matrix (only defined for multiple contrasts)
	Tlabels <- NULL
	# Progressive transformation of Tm, factor by factor
	for (fname in factor.names){
		# f is the current factor vector
		f <- frame[[fname]]
		# n is the factor vector length
		n <- m*nc
		# Remove the corresponding column from the model data frame
		frame[[fname]] <- NULL
		# nc is the number of contrasts for the current factor (updated)
		nc <- ncol(levels[[fname]])
		if (nc > 1L){
			# If there are multiple contrasts, the rows of the model data frame are
			# replicated (by the number of contrasts), and a new column is added to
			# assign a contrast to each group of copied rows. 
			frame[[ncol(frame)+1]] <- 1
			frame <- frame[rep(1:n,nc),]
			frame[[ncol(frame)]] <- rep(1:nc,each=n)
			# Moreover, we create labels to identify what contrast is represented
			# by each row of the final transformation matrix
			new.Tlabels <- paste(fname,as.character(1:nc),sep="")
			if (is.null(Tlabels)){
				Tlabels <- new.Tlabels
			}else{
				# (In case there is more than one factor with multiple contrasts)
				Tlabels <- paste(rep(Tlabels,nc), rep(new.Tlabels,each=length(Tlabels)), sep=":")
			}
		}
		# The same routine as for the first Tm: vector with combined values
		# of the (transformed) model data frame...
		frame_1 <- apply(frame,1,paste,collapse="")
		# ... subset to unique rows
		frame <- unique(frame)
		frame_1.m <- unique(frame_1)
		m <- nrow(frame)/nc
		# And create transformation matrix, depending on the combination of factor
		# levels defined in levels
		kf <- t(levels[[fname]])
		kf <- matrix(rep(kf[,f],each=m),ncol=n)
		To <- matrix(rep(matrix(frame_1,ncol=n,byrow=TRUE),each=m),ncol=n)
		Tu <- matrix(rep(frame_1.m,n),ncol=n)
		Tm <- ((To==Tu) * kf) %*% Tm
	}
	# When the loop is finished, assign labels to final Tm, and return it
	rownames(Tm) <- Tlabels
	return(Tm)
}

## END OF AUXILIARY FUNCTIONS
## MAIN PROCEDURE

## MLM METHOD
testFactors.mlm <- function(model,levels,covariates,offset,terms.formula=~1,inherit.contrasts=FALSE,default.contrasts=c("contr.sum","contr.poly"),idata,icontrasts=default.contrasts,lht=TRUE,...){	
	
	# 1. Make complete list of variables, and between/within-subjects data frames
	if (missing(idata)){
		model.variables <- getPredictors(model)
		within.frame <- NULL
	}else{
		model.variables <- c(getPredictors(model),names(idata))
		# Check for duplicates
		if (any(duplicated(model.variables))) stop("There are redundant variables in the terms of the model")
		within.frame <- idata
	}
	# Include factors (in model$xlevels), with appropriate contrasts
	between.frame <- expand.grid(model$xlevels)
	between.factors <- names(between.frame)
	for (bf in between.factors){
		contrasts(between.frame[[bf]]) <- model$contrasts[[bf]]
	}
	# Define contrasts of factors for testing (not necessarily the same as between.frame)
	# Copy the contrasts explicitly defined in the model frame
	factor.contrasts <- lapply(model.frame(model)[between.factors],"attr","contrasts")
	# Set a value for those that were undefined
	undefined.contrasts <- sapply(factor.contrasts,"is.null")
	if (length(factor.contrasts)>0){
		if (inherit.contrasts){
			factor.contrasts[undefined.contrasts] <- model$contrasts[undefined.contrasts]
		}else{
			are.ordered <- sapply(between.frame,"is.ordered")
			factor.contrasts[undefined.contrasts & are.ordered] <- default.contrasts[2]
			factor.contrasts[undefined.contrasts & !are.ordered] <- default.contrasts[1]
		}
	}

	# Adjust contrasts for the within-subjects frame as well
	within.factors <- names(within.frame)
	for (wf in within.factors){
		if (is.null(attr(within.frame[[wf]], "contrasts"))){
			contrasts(within.frame[[wf]]) <- if (is.ordered(within.frame[[wf]])) icontrasts[2] else icontrasts[1]
		}
	}
	factor.contrasts <- c(factor.contrasts,lapply(within.frame,"attr","contrasts"))
	# Contrast matrices for testing terms
	contrast.mat <- factor.contrasts
	not.mat <- !sapply(contrast.mat,"is.numeric")
	factor.nlevels <- c(lapply(between.frame,"nlevels"),lapply(within.frame,"nlevels"))
	factor.nlevels <- lapply(factor.nlevels,"as.list")
	contrast.mat[not.mat] <- mapply(do.call,factor.contrasts[not.mat],factor.nlevels[not.mat],SIMPLIFY=FALSE)

	# 2. Set the values of numeric predictors (covariates and offset)
	# They are the variables from the model, excluding the response and factors
	numeric.predictors <- model.variables[!(model.variables %in% c(between.factors,within.factors))]
	nv <- setNumericPredictors(model, numeric.predictors,
		if(!missing(covariates)) covariates, if(!missing(offset)) offset)
	covariates <- nv$covariates
	offsetvalue <- nv$offsetvalue

	# 3. Redefine levels into a list of numeric matrices
	# Ignore elements in levels that are not factors of the models
	if (missing(levels)) levels <- NULL	else{
		# See if contrastCoefficients must be applied
		if (is.null(names(levels))){
			levels <- do.call("contrastCoefficients", c(levels, list(data=model$model)))
		}else{
			gotnames <- as.logical(sapply(names(levels), nchar))
			levels[!gotnames] <- do.call("contrastCoefficients", c(levels[!gotnames], list(data=model$model)))
		}
		valid.levels <- names(levels) %in% c(between.factors, within.factors)
		if (!all(valid.levels)) warning("Some variables in levels are not model factors, and will be ignored.")
		levels <- levels[valid.levels]
	}
	# Search factors of levels in between.frame and within.frame, and transform the level
	# combinations into a numeric matrix format
	interacting.factors <- names(levels)
	between.factors <- between.factors[between.factors %in% interacting.factors]
	levels <- checkFactors(between.frame,between.factors,levels)
	within.factors <- within.factors[within.factors %in% interacting.factors]
	levels <- checkFactors(within.frame,within.factors,levels)

	# 4. Analyse each term of the terms.formula
	# Stop if the variables of the formula are not in the model
	user.variables <- all.vars(terms.formula)
	if (!all(user.variables %in% model.variables)) stop("The variables in terms.formula are not coherent with the model")
	# Split the formula in terms, and terms in variables, and add Intercept if suitable
	user.terms <- terms(terms.formula)
	user.terms.labels <- labels(user.terms)
	split.user.terms <- strsplit(user.terms.labels,":")
	if (attr(user.terms,"intercept")){
		user.terms.labels <- c("(Intercept)",user.terms.labels)
		split.user.terms <- c("(Intercept)",split.user.terms)
	}
	# term is a list with the variables that will change in each analysis
	# test.result will contain the results of the analysis 
	term <- vector("list",7L)
	names(term) <- c("num.vars","fac.vars","levels","between.factors","within.factors","covariates","offset")
	test.result <- vector("list", length(split.user.terms))
	names(test.result) <- user.terms.labels
	# Repeat for each term:
	for (i in seq(length(split.user.terms))){
		# Copy default values of term
		term$num.vars <- term.vars <- split.user.terms[[i]]
		term$fac.vars <- character()
		term$levels <- levels
		term$between.factors <- between.factors
		term$within.factors <- within.factors
		term$covariates <- covariates
		term$offset <- offsetvalue
		# Take factors in term, add their contrast matrices to the copy of levels,
		# and then remove them from the list of term variables,
		# which now will only contain the covariates of the term;
		# first for between-subjects predictors
		replace.factors <- term.vars[term.vars %in% names(between.frame)]
		if (length(replace.factors) > 0){
			term$levels[replace.factors] <- contrast.mat[replace.factors]
			term$between.factors <- unique(c(between.factors,replace.factors))
			term$num.vars <- term$num.vars[-match(replace.factors,term$num.vars)]
			term$fac.vars <- replace.factors
		}
		# and the same for within-subjects
		replace.factors <- term.vars[term.vars %in% names(within.frame)]
		if (length(replace.factors) > 0){
			term$levels[replace.factors] <- contrast.mat[replace.factors]
			term$within.factors <- unique(c(within.factors,replace.factors))
			term$num.vars <- term$num.vars[-match(replace.factors,term$num.vars)]
			term$fac.vars <- c(term$fac.vars,replace.factors)
		}
		if (length(term$num.vars)==0) term$num.vars <- "(Intercept)"
		# Take covariates in term, and change their value to 1 in the argument covariates
		if (term$num.vars[1] != "(Intercept)") term$covariates[term$num.vars] <- 1
		# Run basic analysis
		test.result[[user.terms.labels[i]]] <- testFactorsOnTerm.mlm(model,term,numeric.predictors,between.frame,within.frame,lht,...)
	}
	# Check if the tested contrasts of factors are orthogonal to the intercept
	factor.contrasts <- factor.contrasts[names(factor.contrasts) %in% user.variables]
	contrast.mat <- contrast.mat[names(contrast.mat) %in% user.variables]
	not.orthogonal <- sapply(contrast.mat,function(x) any(as.logical(round(colSums(x),digits=getOption("digits")))))
	if (any(not.orthogonal)) warning(paste("Contrasts are not orthogonal for factor(s):",paste(names(contrast.mat)[not.orthogonal],collapse=", ")))
	result <- list(call=match.call(,sys.call(1L)),model.call=model$call,levels=levels,factor.contrasts=factor.contrasts,covariates=covariates,terms=test.result)
	class(result) <- c("testFactors.mlm","testFactors")
	return(result)
}

testFactorsOnTerm.mlm <- function(model,term,numeric.predictors,between.frame,within.frame,lht,...){
	# 1. Redefine model data frame
	# Include columns with averages or user-defined values for the numeric predictors
	between.frame <- addPredictorsToFrame(between.frame, numeric.predictors, term)
	
	# 2. Preliminary definition of the Linear Hypothesis matrix (L)
	L <- defineLHT(model,term,between.frame)
	# (add offset column if available)
	if (term$offset != 0) L <- cbind(L, term$offset)
	
	# 3. Transformed Linear Hypothesis (L) and response transformation (P) matrices
	if (length(term$between.factors > 0)){
		L <- makeT(between.frame,term$between.factors,term$levels) %*% L
	}else{
		L <- t(as.matrix(colSums(L)/nrow(L)))
	}
	if (length(term$within.factors) > 0){
		P <- t(makeT(within.frame,term$within.factors,term$levels))
	}else{
		P <- if (is.null(term$within.factors)) NULL else matrix(rep(1/nrow(within.frame),nrow(within.frame)))
	}
	# (and separate offset)
	if (term$offset){
		offset_effect <- L[,ncol(L)]
		L <- L[,-ncol(L), drop=FALSE]
	}else offset_effect <- 0
	
	# 4. Result, consisting in:
	#   levels: numeric matrix values of levels
	#   adjusted.values: table of adjusted means for the tested interactions
	#   covmat: covariance matrix the adjusted values
	#   std.error: standard error of adjusted values (from covmat)
	#   test: test value, from LinearHypothesis
	adjusted.values <- L %*% model$coefficients + offset_effect
	Ln <- diag(ncol(model$coefficients)) %x% L
	covmat <- Ln %*% vcov(model) %*% t(Ln)
	if (!is.null(P)){
		adjusted.values <- adjusted.values %*% P
		Pm <- t(P) %x% diag(nrow(L))
		covmat <- Pm %*% covmat %*% t(Pm)
		rownames(P) <- colnames(model$coefficients)
	}
	# define dimnames(covmat)?
	std.error <- matrix(sqrt(diag(covmat)),ncol=ncol(adjusted.values))
	dimnames(std.error) <- dimnames(adjusted.values)
	result <- list(numeric.variables=paste(term$num.vars,sep=":"),factor.variables=term$fac.vars,hypothesis.matrix=L,P=P,adjusted.values=adjusted.values,covmat=covmat,std.error=std.error)
	if (lht) result <- c(result,list(test=try(linearHypothesis(model,L,P=P,...),silent=TRUE)))
	return(result)
}

## DEFAULT METHOD
testFactors.default <- function(model,levels,covariates,offset,terms.formula=~1,inherit.contrasts=FALSE,default.contrasts=c("contr.sum","contr.poly"),lht=TRUE,...){	

	# 1. Make complete list of predictor variables, and extract factors
	predictors <- getPredictors(model)
	X <- getModelFrame(model)[predictors]
	are.factors <- as.logical(sapply(X, is.factor))
	factor.names <- predictors[are.factors] 
	# Factor data frame, with appropriate contrasts
	factor.frame <- expand.grid(lapply(X[are.factors],"levels"))
	for (f in factor.names){
		contrasts(factor.frame[[f]]) <- getContrasts(model)[[f]]
	}
	# Define contrasts of factors for testing (not necessarily the same as factor.frame)
	# Copy the contrasts explicitly defined in the model frame
	factor.contrasts <- lapply(X[are.factors],"attr","contrasts")
	# Set a value for those that were undefined
	undefined.contrasts <- sapply(factor.contrasts,"is.null")
	if (length(factor.contrasts)>0){
		if (inherit.contrasts){
			factor.contrasts[undefined.contrasts] <- getContrasts(model)[undefined.contrasts]
		}else{
			are.ordered <- sapply(factor.frame,"is.ordered")
			factor.contrasts[undefined.contrasts & are.ordered] <- default.contrasts[2]
			factor.contrasts[undefined.contrasts & !are.ordered] <- default.contrasts[1]
		}
	}
	# Contrast matrices for testing terms
	contrast.mat <- factor.contrasts
	not.mat <- !sapply(contrast.mat,"is.numeric")
	factor.nlevels <- lapply(lapply(factor.frame,"nlevels"),"as.list")
	contrast.mat[not.mat] <- mapply(do.call,factor.contrasts[not.mat],factor.nlevels[not.mat],SIMPLIFY=FALSE)
	
	# 2. Set the values of numeric predictors (covariates and offset)
	numeric.predictors <- predictors[!are.factors]
	nv <- setNumericPredictors(model, numeric.predictors,
		if(!missing(covariates)) covariates, if(!missing(offset)) offset)
	covariates <- nv$covariates
	offsetvalue <- nv$offsetvalue
	
	# 3. Redefine levels into a list of numeric matrices
	# Ignore elements in levels that are not factors of the models
	if (missing(levels)) levels <- NULL else{
		# See if contrastCoefficients must be applied
		if (is.null(names(levels))){
			levels <- do.call("contrastCoefficients", c(levels, list(data=getModelFrame(model))))
		}else{
			gotnames <- as.logical(sapply(names(levels), nchar))
			levels[!gotnames] <- do.call("contrastCoefficients", c(levels[!gotnames], list(data=getModelFrame(model))))
		}
		valid.levels <- names(levels) %in% factor.names
		if (!all(valid.levels)) warning("Some variables in levels are not model factors, and will be ignored.")
		levels <- levels[valid.levels]
	}
	# Search factors of levels in factor.frame, and transform the level
	# combinations into a numeric matrix format
	interacting.factors <- names(levels)
	factor.names <- factor.names[factor.names %in% interacting.factors]
	levels <- checkFactors(factor.frame,factor.names,levels)
	
	# 4. Analyse each term of the terms.formula
	# Stop if the variables of the formula are not in the model
	user.variables <- as.character(attr(terms(terms.formula),"variables")[-1])
	if (!all(user.variables %in% predictors)) stop("The variables in terms.formula are not coherent with the model")
	# Split the formula in terms, and terms in variables, and add Intercept if suitable
	user.terms <- terms(terms.formula)
	user.terms.labels <- labels(user.terms)
	split.user.terms <- strsplit(user.terms.labels,":")
	if (attr(user.terms,"intercept")){
		user.terms.labels <- c("(Intercept)",user.terms.labels)
		split.user.terms <- c("(Intercept)",split.user.terms)
	}
	# term is a list with the variables that will change in each analysis
	# test.result will contain the results of the analysis 
	term <- vector("list",6L)
	names(term) <- c("num.vars","fac.vars","levels","factor.names","covariates","offset")
	test.result <- vector("list", length(split.user.terms))
	names(test.result) <- user.terms.labels
	# Repeat for each term:
	for (i in 1:length(split.user.terms)){
		# Copy default values of term
		term$num.vars <- term.vars <- split.user.terms[[i]]
		term$fac.vars <- character()
		term$levels <- levels
		term$factor.names <- factor.names
		term$covariates <- covariates
		term$offset <- offsetvalue
		# Take factors in term, add their contrast matrices to the copy of levels,
		# and then remove them from the list of term variables,
		# which now will only contain the covariates of the term;
		replace.factors <- term.vars[term.vars %in% names(factor.frame)]
		if (length(replace.factors) > 0){
			term$levels[replace.factors] <- contrast.mat[replace.factors]
			term$factor.names <- unique(c(factor.names,replace.factors))
			term$num.vars <- term$num.vars[-match(replace.factors,term$num.vars)]
			term$fac.vars <- replace.factors
		}
		if (length(term$num.vars)==0) term$num.vars <- "(Intercept)"
		# Take covariates in term, and change their value to 1 in the argument covariates
		if (term$num.vars[1] != "(Intercept)") term$covariates[term$num.vars] <- 1
		# Run basic analysis
		test.result[[user.terms.labels[i]]] <- testFactorsOnTerm.default(model,term,numeric.predictors,factor.frame,lht,...)
	}
	# Check if the tested contrasts of factors are orthogonal to the intercept
	factor.contrasts <- factor.contrasts[names(factor.contrasts) %in% user.variables]
	contrast.mat <- contrast.mat[names(contrast.mat) %in% user.variables]
	not.orthogonal <- sapply(contrast.mat,function(x) any(as.logical(round(colSums(x),digits=getOption("digits")))))
	if (any(not.orthogonal)) warning(paste("Contrasts are not orthogonal for factor(s):",paste(factor.names[not.orthogonal],collapse=",")))
	result <- list(call=match.call(,sys.call(1L)),model.call=getCall(model),levels=levels,factor.contrasts=factor.contrasts,covariates=covariates,terms=test.result)
	class(result) <- "testFactors"
	return(result)
}

testFactorsOnTerm.default <- function(model,term,numeric.predictors,factor.frame,lht,...){
	# 1. Redefine model data frame
	# Include columns with averages or user-defined values for the numeric predictors
	factor.frame <- addPredictorsToFrame(factor.frame, numeric.predictors, term)
	
	# 2. Preliminary definition of the Linear Hypothesis matrix (L)
	L <- defineLHT(model,term,factor.frame)
	# (add offset column if available)
	if (term$offset != 0) L <- cbind(L, term$offset)

	# 3. Transformed Linear Hypothesis (L)
	if (length(term$factor.names > 0)){
		L <- makeT(factor.frame,term$factor.names,term$levels) %*% L
	}else{
		L <- t(as.matrix(colSums(L)/nrow(L)))
	}
	# (and separate offset)
	if (term$offset != 0){
		offset_effect <- L[,ncol(L)]
		L <- L[,-ncol(L), drop=FALSE]
	}else offset_effect <- 0

	# 4. Result, consisting in:
	#   levels: numeric matrix values of levels
	#   adjusted.values: table of adjusted means for the tested interactions
	#   covmat: covariance matrix the adjusted values
	#   std.error: standard error of adjusted values (from covmat)
	#   test: test value, from LinearHypothesis
	adjusted.values <- L %*% getCoef(model) + offset_effect
	covmat <- L %*% vcov(model) %*% t(L)
	# define dimnames(covmat)?
	std.error <- matrix(sqrt(diag(covmat)), ncol=ncol(adjusted.values))
	dimnames(std.error) <- dimnames(adjusted.values)
	result <- list(numeric.variables=paste(term$num.vars,sep=":"),factor.variables=term$fac.vars,hypothesis.matrix=L,adjusted.values=adjusted.values,covmat=covmat,std.error=std.error)
	if (lht) result <- c(result,list(test=try(linearHypothesis(model,L,...),silent=TRUE)))
	return(result)
}

## LM METHOD (equal to default, but with a result of a different class, just for summary purposes)
testFactors.lm <- function(model,...){
	result <- testFactors.default(model,...)
	class(result) <- c("testFactors.lm","testFactors")
	return(result)
}

## GLM METHOD (equal to default, but adjusted means are transformed)
testFactors.glm <- function(model, ..., link=FALSE){
	result <- testFactors.default(model,...)
	attr(result,"means") <- if (link) "link" else "mean"
	if (!link){
		terms.with.means <- (lapply(result$terms,"[[","numeric.variables")=="(Intercept)")
		for (term in which(terms.with.means)) result$terms[[term]]$adjusted.values <- getFamily(model)$linkinv(result$terms[[term]]$adjusted.values)
	}
	class(result) <- c("testFactors.glm","testFactors")
	return(result)
}

## LME METHOD (equal to default, but with a result of a different class, just for summary purposes)
testFactors.lme <- function(model, ...){
	result <- testFactors.default(model,...)
	class(result) <- c("testFactors.lme","testFactors")
	return(result)
}

## MER METHOD (equal to default, but adjusted means are transformed if suitable)
testFactors.mer <- function(model, ..., link=FALSE){
	if (!as.logical(model@dims["LMM"]) && length(model@muEta) == 0){
	stop("Nonlinear Mixed Models are not supported.")
	}
	result <- testFactors.default(model,...)
	if (length(model@muEta > 0)){
		attr(result,"means") <- if (link) "link" else "mean"
		if (!link){
			terms.with.means <- (lapply(result$terms,"[[","numeric.variables")=="(Intercept)")
			for (term in which(terms.with.means)){
				result$terms[[term]]$adjusted.values <- getFamily(model)$linkinv(result$terms[[term]]$adjusted.values)
			}
		}
	}
	class(result) <- c("testFactors.mer","testFactors")
	return(result)
}

## MERMOD METHOD (equal to default, but adjusted means are transformed if suitable)
testFactors.merMod <- function(model, ..., link=FALSE){
	if (is(model, "nlmerMod")){stop("Nonlinear Mixed Models are not supported.")}
	result <- testFactors.default(model,...)
	if (is(model, "lmerMod")) attr(result,"means") <- "mean"
	if (is(model, "glmerMod")){
		attr(result,"means") <- if (link) "link" else "mean"
		if (!link){
			terms.with.means <- (lapply(result$terms,"[[","numeric.variables")=="(Intercept)")
			for (term in which(terms.with.means)){
				result$terms[[term]]$adjusted.values <- getFamily(model)$linkinv(result$terms[[term]]$adjusted.values)
			}
		}
	}
	class(result) <- c("testFactors.merMod","testFactors")
	return(result)
}


testFactors <- function(model,...){UseMethod("testFactors")}

## PRINT AND SUMMARY METHODS

print.testFactors <- function(x,digits=getOption("digits"),...){
	cat("\nCall:",deparse(x$call),"\n")
	mean.label <- if (as.character(x$model.call)[1] %in% c("glm", "glmer") && attr(x,"means")=="link") "link function" else "mean"
	for (i in seq(length(x$terms))){
		cat("\nTerm",names(x$terms)[i],"\n")
		term <- x$terms[[i]]
		if (term$numeric.variables[1] == "(Intercept)") cat("\nAdjusted",mean.label) else cat("\nAdjusted slope for",paste(term$numeric.variables,collapse=":"))
		if (length(term$factor.variables)>0) cat(" at contrasts of",paste(term$factor.variables,collapse=", "))
		cat(":\n")
		if (dim(term$adjusted.values)[2]==1) dimnames(term$adjusted.values)[2] <- list("")
		print(drop(term$adjusted.values),digits=digits,...)
		cat("\nStd. Error")
		if (as.character(x$model.call)[1] %in% c("glm", "glmer") && term$numeric.variables[1] == "(Intercept)"){
			cat(" of link function")
		}
		cat(":\n")
		if (dim(term$std.error)[2]==1) dimnames(term$std.error)[2] <- list("")
		print(drop(term$std.error),digits=digits,...)
		if ("test" %in% names(term) && class(term$test)[1] != "try-error"){
			cat("\n") #cat("\nLinear hypothesis test:\n")
			print(term$test,digits=digits,...)
		}
		cat("------\n")
	}
	invisible(x)
}

summary.testFactors <- function(object,predictors=TRUE,matrices=TRUE,covmat=FALSE,...){
	sobject <- list()
	sobject$model.call <- object$model.call
	attr(sobject,"means") <- attr(object,"means")
	# Define values of predictors, if requested
	if (predictors){
		sobject$levels <- object$levels
		sobject$factor.contrasts <- object$factor.contrasts
		sobject$covariates <- object$covariates
	}
	# Define test details (matrices L and P), if requested
	if (matrices){
		# Matrices L, bound by rows
		lh.matrices <- lapply(object$terms,"[[","hypothesis.matrix")
		lh.labels1 <- rep(names(lh.matrices), sapply(lh.matrices,"nrow"))
		lh.matrices <- do.call("rbind",lh.matrices)
		# If the original matrices have row names (there are various contrasts), copy those names in the labels
		lh.labels2 <- rownames(lh.matrices)
		rownames(lh.matrices) <- NULL
		if (length(lh.labels2)>0){
			# In secondary labels, append "|" to separate from the primary label
			lh.labels2[nchar(lh.labels2)>0] <- paste("|",lh.labels2[nchar(lh.labels2)>0])
			sobject$hypothesis.matrix <- as.data.frame(lh.matrices,row.names=paste(lh.labels1,lh.labels2))
		}else{
			sobject$hypothesis.matrix <- as.data.frame(lh.matrices,row.names=lh.labels1)
		}
		# Do the same for matrice P (if they exist), bound by columns and then transposed
		p.matrices <- lapply(object$terms,"[[","P")
		if (!sapply(p.matrices,"is.null")[1]){
			p.labels1 <- rep(names(p.matrices), sapply(p.matrices,"ncol"))
			p.matrices <- do.call("cbind",p.matrices)
			p.labels2 <- colnames(p.matrices)
			colnames(p.matrices) <- NULL
			if (length(p.matrices)>0){
				if (length(p.labels2)>0){
					p.labels2[nchar(p.labels2)>0] <- paste("|",p.labels2[nchar(p.labels2)>0])
					sobject$P <- as.data.frame(t(p.matrices),row.names=paste(p.labels1,p.labels2))
				}else{
					sobject$P <- as.data.frame(t(p.matrices),row.names=p.labels1)
				}
			}
		}
	}
	# Build list of adjusted values and ANOVA table
	sobject$adjusted.values <- lapply(object$terms,"[[","adjusted.values")
	sobject$std.error <- lapply(object$terms,"[[","std.error")
	if (covmat) sobject$covmat <- lapply(object$terms,"[[","covmat")
	for (term.label in names(object$terms)){
		# The numeric and factor variables referred by the adjusted values are assigned to attributes
		attr(sobject$adjusted.values[[term.label]],"numeric.variables") <- object$terms[[term.label]]$numeric.variables
		attr(sobject$adjusted.values[[term.label]],"factor.variables") <- object$terms[[term.label]]$factor.variables
	}
	sobject$anova.table <- anova(object)
	class(sobject) <- "summary.testFactors"
	return(sobject)
}

anova.testFactors <- function(object,...){
	tests <- lapply(object$terms,"[[","test")
	# The ANOVA table will be built for terms where linearHypothesis was successfully used
	successful.tests <- sapply(tests,function(x) !is.null(x) && class(x)[1]!="try-error")
	anova.table <- matrix(NA,sum(successful.tests),3)
	rownames(anova.table) <- names(object$terms[successful.tests])
	for (term.label in names(object$terms)){
		# Create ANOVA table, copying values (Df, Chisq/F statistic, and P) from the test result
		if (successful.tests[term.label]){
			lht.result <- object$terms[[term.label]]$test
			anova.table[term.label,1:3] <- as.matrix(lht.result[2, seq(to=ncol(lht.result),length=3)])
		}
	}
	if (nrow(anova.table) > 0){
		# Add a row with Df of the residual if the information is available
		if (ncol(lht.result) > 3) {anova.table <- rbind(anova.table, Residual=c(lht.result[[2,1]],NA,NA))}
		colnames(anova.table) <- colnames(lht.result)[seq(to=ncol(lht.result),length=3)]
		anova.table <- structure(as.data.frame(anova.table),
			heading = paste(colnames(anova.table)[ncol(anova.table)-1], " Test: ", sep = ""), class = c("anova", "data.frame"))
	}
	return(anova.table)
}

anova.testFactors.lm <- function(object,predictors=TRUE,matrices=TRUE,...){
	tests <- lapply(object$terms,"[[","test")
	# The ANOVA table will be built for terms where linearHypothesis was successfully used
	successful.tests <- sapply(tests,function(x) !is.null(x) && class(x)[1]!="try-error")
	anova.table <- matrix(NA,sum(successful.tests),4)
	rownames(anova.table) <- names(object$terms[successful.tests])
	for (term.label in names(object$terms)){
		# Create ANOVA table, copying values (Df, SS, F/Chisq statistic, and P) from the test result
		if (successful.tests[term.label]){
			lht.result <- object$terms[[term.label]]$test
			anova.table[term.label,1:4] <- as.matrix(lht.result[2,3:6])
		}
	}
	if (nrow(anova.table) > 0){
		# Add a row with Df of the residual
		anova.table <- rbind(anova.table,
			Residual=c(as.matrix(lht.result[2,1:2]),NA,NA))
		colnames(anova.table) <- colnames(lht.result)[3:6]
		anova.table <- structure(as.data.frame(anova.table), heading = paste(colnames(anova.table)[3], " Test: ", sep = ""), class = c("anova", "data.frame"))
	}
	return(anova.table)
}

anova.testFactors.mlm <- function(object,...){
	tests <- lapply(object$terms,"[[","test")
	# The ANOVA table will be built for terms where linearHypothesis was successfully used
	successful.tests <- sapply(tests,function(x) !is.null(x) && class(x)[1]!="try-error")
	anova.table <- matrix(NA,sum(successful.tests),4)
	rownames(anova.table) <- names(object$terms[successful.tests])
	for (term.label in names(object$terms)){
		# Create ANOVA table, copied from print.linearHypothesis.mlm
		if (successful.tests[term.label]){
			lht.result <- object$terms[[term.label]]$test
			test <- lht.result$test[1]
			SSPE.qr <- qr(lht.result$SSPE)
			eigs <- Re(eigen(qr.coef(SSPE.qr, lht.result$SSPH), symmetric = FALSE)$values)
			anova.table[term.label,1:4] <- switch(test,
				"Pillai" = stats_Pillai(eigs, lht.result$df, lht.result$df.residual),
				"Wilks" = stats_Wilks(eigs, lht.result$df, lht.result$df.residual),
				"Hotelling-Layley" = stats_HL(eigs, lht.result$df, lht.result$df.residual),
				"Roy" = stats_Roy(eigs, lht.result$df, lht.result$df.residual))
		}
	}
	# Complete ANOVA table, like in print.linearHypothesis.mlm
	if (nrow(anova.table) > 0){
		ok <- anova.table[, 2] >= 0 & anova.table[, 3] > 0 & anova.table[, 4] > 0
		ok <- !is.na(ok) & ok
		anova.table <- cbind(lht.result$df, anova.table, pf(anova.table[ok, 2], anova.table[ok, 3], anova.table[ok, 4], lower.tail = FALSE))
		colnames(anova.table) <- c("Df", "test stat", "approx F", "num Df", 
		"den Df", "Pr(>F)")
		anova.table <- structure(as.data.frame(anova.table), heading = paste("Multivariate Test", if (nrow(anova.table) > 1) "s", ": ", test, " test statistic", sep = ""), class = c("anova", "data.frame"))
	}
	return(anova.table)
}

print.summary.testFactors <- function(x,digits=getOption("digits"),...){
	cat("\nAdjusted values for factor combinations in model:\n",deparse(x$model.call),"\n")
	elements <- names(x)
	# List of levels and the other details of predictors, if had been requested by summary
	if ("levels" %in% elements){
		header <- "\nValues of predictor variables.\n"
		if (length(x$levels)>0){
			cat(header,"\nSpecified combinations of factor levels:\n")
			header <- "\n---\n"
			for (nfac in names(x$levels)){
				cat("\n") #cat("\n",nfac,"\n")
				fac <- x$levels[[nfac]]
				print(t(fac),digits=digits,...)
			}
		}
		if (length(x$factor.contrasts)>0){
			cat(header,"\nDefault list of factor contrasts:\n")
			header <- "\n---\n"
			factor.contrasts <- sapply(x$factor.contrasts,"[")
			if (class(factor.contrasts)=="character"){
				print(factor.contrasts,quote=FALSE)
			}else{
				for (nfac in names(factor.contrasts)){
					cat("\n",nfac,"\n")
					fac <- factor.contrasts[[nfac]]
					if (class(fac)=="character"){
						cat(":",fac,"\n")
					}else{
						print(fac)
						cat("\n")
					}
				}
			}
		}
		if (length(x$covariates)){
			cat(header,"\nSpecified values of covariates:\n")
			print(x$covariates)
		}
	}
	# Matrices L and P, if requested (and exist)
	if ("hypothesis.matrix" %in% elements){
		cat("------\n\nLinear hypothesis matrix\n\n")
		print(x$hypothesis.matrix)
		cat("\n")
	}
	if (("P" %in% elements) && (length(x$P)>0)){
		cat("\n---\nResponse transformation matrix\n\n")
		print(x$P)
		cat("\n")
	}
	# Adjusted values and standard errors (plus covariance matrices, if requested and exist)
	cat("------\n\nAdjusted values\n")
	mean.label <- if (as.character(x$model.call)[1] %in% c("glm", "glmer") && attr(x,"means")=="link") "link function" else "mean"
	for (n in names(x$adjusted.values)){
		cat("\nTerm",n,"\n")
		mat <- x$adjusted.values[[n]]
		term_name <- attr(mat,"numeric.variables")[1]
		if (term_name == "(Intercept)") cat("\nAdjusted",mean.label) else cat("\nAdjusted slope for",paste(attr(mat,"numeric.variables"),collapse=":"))
		if (length(attr(mat,"factor.variables"))>0) cat(" at contrasts of",paste(attr(mat,"factor.variables"),collapse=", "))
		cat(":\n")
		attr(mat,"numeric.variables")<-NULL
		attr(mat,"factor.variables")<-NULL
		if (dim(mat)[2]==1) dimnames(mat)[2] <- list("")
		print(drop(mat),digits=digits,...)
		cat("\nStandard error")
		if (as.character(x$model.call)[1] %in% c("glm", "glmer") && term_name == "(Intercept)"){
			cat(" of link function")
		}
		cat(":\n")
		mat <- x$std.error[[n]]
		if (dim(mat)[2]==1) dimnames(mat)[2] <- list("")
		print(drop(mat),digits=digits,...)
		if ("covmat" %in% elements){
			cat("\nVariance-covariance matrix:\n")
			mat <- x$covmat[[n]]
			if (dim(mat)[2]==1) dimnames(mat)[2] <- list("")
			print(drop(mat),digits=digits,...)
		}
		cat("---\n")
	}
	# ANOVA table (if exists)
	if (nrow(x$anova.table)>0){
		cat("\n------\n") #cat("---\n\nLinear hypothesis tests\n\n")
		print(x$anova.table)
	}
	cat("------------\n")
	invisible(x)
}
