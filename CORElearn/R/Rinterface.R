versionCore <- function()
{
	tmp <- .C("versionCore",
			libVersion = paste(rep(".",times=256),collapse=""),
			PACKAGE="CORElearn"
	)
	tmp$libVersion
}
destroyModels <-function(model=NULL)
{
	if (is.null(model)) {
		destroyCore()
		initCore()
	}
	else {
		for (i in seq(along=model)) {
			modelID <- model[i]$modelID
			tmp <- .C("destroyOneCoreModel", as.integer(modelID), PACKAGE="CORElearn" )
		}
	}
	invisible(NULL) 
}
CoreModel <- function(formula, data, model=c("rf","rfNear","tree","knn","knnKernel","bayes","regTree"), ..., costMatrix=NULL)
{
	# check formula or response index or reponse name
	if (inherits(formula,"formula")) {
		dat <- model.frame(formula, data=data, na.action=na.pass)
		trms <- attr(dat,"terms")
		attributes(trms) <- NULL
		formulaExpanded <- as.formula(trms)
	} 
	else {
		if (is.numeric(formula)) {
			if (formula == round(formula)) {# index of response variable
				classIdx <- formula
				className <- names(data)[classIdx]
			}
			else  stop("The first argument must be a formula or prediction column name or prediction column index.")
		}
		else if (is.character(formula)) { # name of response variable
			classIdx <- match(formula, names(data))
			if (length(classIdx) != 1 || is.na(classIdx)) 
				stop("The first argument must be a formula or prediction column name or prediction column index.")
			className <- names(data[classIdx])
		}
		else stop("The first argument must be a formula or prediction column name or prediction column index.")
		
		dat <- data.frame(data[, classIdx], data[, -classIdx, drop=FALSE])
		names(dat)[1] <- className
		# get formula explicitly to allow storage of all terms and their manipulation
		frml <- paste(className, "~",paste(names(dat)[-1], sep="+",collapse="+"),sep="") 
		formulaExpanded <- as.formula(frml)   
	}
	
	model <- match.arg(model) ; 
	isRegression <- model == "regTree"
	if (isRegression && !is.null(costMatrix))
		warning("For regression problems parameter costMatrix is ignored.");  
	
	if (!isRegression && !inherits(dat[[1]],"factor")) {
		dat[[1]] <- factor(dat[[1]]);
		cat("Changing dependent variable to factor with levels:",levels(dat[[1]]),"\n");
		warning("Possibly this is an error caused by regression formula and classification model or vice versa.")
	}
	class.lev <- levels(dat[[1]]);
	noClasses <- length(class.lev);
	if (!isRegression && is.null(costMatrix)) {
		## create and fill uniform costs matrix
		costMatrix <- 1 - diag(noClasses);
	}
	
	aux <- prepare.Data(dat, formulaExpanded, dependent=TRUE,numericAsOrdered=FALSE,skipNAcolumn=TRUE,skipEqualColumn=TRUE);
	discnumvalues <- aux$discnumvalues;
	discdata <- aux$discdata;
	numdata <- aux$numdata;
	discAttrNames <- dimnames(discdata)[[2]]
	discValCompressed <- aux$disccharvalues
	discValues <- aux$discValues
	numAttrNames <- dimnames(numdata)[[2]]
	skipNames <- dimnames(aux$skipmap)
	discmap <- aux$discmap;
	nummap <- aux$nummap;    
	skipmap <- aux$skipmap
	options <- prepare.Options(...);
	checkOptionsValues(options) ;
	optRemain <- checkModelOptions(model, options) ;
	if (length(optRemain) > 0) warning("Unused options:", paste(names(optRemain), collapse=", "));
	options <- convert.Options(options)
	options <- c(options, action=model)
	tmp <- .C("buildCoreModel",
			noInst = aux$noInst,
			noDiscrete = ncol(discdata),
			noDiscreteValues = as.integer(discnumvalues),
			discreteData = as.integer(discdata), # vector of length noInst*noDiscrete, column-wise
			noNumeric = ncol(numdata),
			numericData = as.double(numdata), # vector of length noInst*noNumeric, column-wise
			costs = as.double(costMatrix),
			discAttrNames = as.character(discAttrNames),
			discValNames = as.character(discValCompressed),
			numAttrNames = as.character(numAttrNames),
			numOptions = length(options),
			optionsName = names(options),
			optionsVal = options,
			modelID = integer(1),
			noClasses = integer(1),
			priorClassProb = numeric(256),
			avgTrainPrediction=numeric(1),
			NAOK=TRUE,
			PACKAGE="CORElearn"
	);
	if (tmp$modelID == -1) {
		return(NULL)
	}
	res <- list(modelID=tmp$modelID, class.lev=class.lev, model=model, formula=aux$formulaOut,
			noClasses = tmp$noClasses, priorClassProb = tmp$priorClassProb[1:tmp$noClasses],
			avgTrainPrediction = tmp$avgTrainPrediction,
			noNumeric = tmp$noNumeric, noDiscrete=tmp$noDiscrete, discAttrNames = discAttrNames,
			discValNames = discValues, numAttrNames = numAttrNames,   
			discmap = discmap, nummap = nummap, skipmap = skipmap
	)
	class(res) <- "CoreModel"
	res
}
predict.CoreModel <- function(object, newdata, ..., costMatrix=NULL,  type=c("both","class","probability"))
{
	model <- object
	#rm(object)
	type<-match.arg(type)
	modelID <- model$modelID;
	isRegression <- model$noClasses == 0
	noClasses <- model$noClasses;
	class.lev <- model$class.lev;
	#terms <- delete.response(model$terms);
	newFormula <- reformulate(all.vars(model$formula)[-1])
	#newdata <- as.data.frame(newdata)
	#dat <- model.frame(model$formula, data=newdata, na.action=na.pass);
	dat <- model.frame(newFormula, data=newdata, na.action=na.pass);
	aux <- prepare.Data(dat, model$formula, dependent=FALSE,class.lev=class.lev, numericAsOrdered=FALSE,skipNAcolumn=FALSE, skipEqualColumn=FALSE);
	#aux <- prepare.Data(dat[-1], model$formula, dependent=FALSE,class.lev, skipNAcolumn=FALSE, skipEqualColumn=FALSE);
	noInst <- aux$noInst
	discnumvalues <- aux$discnumvalues;
	discdata <- aux$discdata;
	numdata <- aux$numdata;
	if (!isRegression && is.null(costMatrix)) 
		costMatrix <- 1 - diag(noClasses); ## create and fill uniform costs matrix
	options <- prepare.Options(...);
	checkOptionsValues(options) ;
	optRemain <- checkPredictOptions(model, options)
	if (length(optRemain) > 0) warning("Unused options:", paste(names(optRemain), collapse=", "));       
	
	tmp <- .C("predictWithCoreModel",
			as.integer(modelID),
			noInst = noInst,
			discreteData = as.integer(discdata), # vector of length noInst*noDiscrete, columnwise
			numericData = as.double(numdata), # vector of length noInst*noNumeric, columnwise
			costs = as.double(costMatrix),
			predicted = integer(noInst),
			prob = double(noInst*noClasses),
			predictedReg = double(noInst),
			numOptions = length(options),
			optionsName = names(options),
			optionsVal = options,
			NAOK=TRUE,
			PACKAGE="CORElearn"
	);
	if (model$model == "regTree") {
		returnList <- tmp$predictedReg
	}
	else {
		code <- tmp$predicted;
		code[code==0] <- NA;
		pred <- factor(class.lev[code],levels=class.lev);
		prob <- matrix(tmp$prob, nrow=noInst, ncol=noClasses,dimnames=list(NULL,class.lev));
		if (type == "both")
			returnList <- list(class=pred,probabilities=prob)
		else if (type=="class")
			returnList <- pred
		else if (type == "probability")
			returnList <- prob
	}
	returnList
}

display<- function(x, format=c("screen","dot")) UseMethod("display", x)

display.CoreModel <- function(x, format=c("screen","dot")) {
	format<-match.arg(format)
	if (x$model %in% c("tree","knn","knnKernel","bayes","regTree")){
		if (format=="screen")
			treeStr <- .Call("printTree2R", as.integer(x$modelID), PACKAGE="CORElearn")
		else if (format == "dot")
			treeStr <- .Call("printTreeDot2R", as.integer(x$modelID), PACKAGE="CORElearn")
	}
	else {
		warning("The model provided is not of appropriate type for this visualization.");
		treeStr <- ""
	}
	cat(treeStr)
	invisible(treeStr)	
} 

plot.CoreModel<-function(x, trainSet, rfGraphType=c("attrEval","outliers","scaling", "prototypes","attrEvalCluster"), clustering=NULL,...) 
{    
	rfGraphType<-match.arg(rfGraphType)
	# regression or decision tree
	if (x$model == "regTree" || x$model == "tree"){           
		rmodel <- getRpartModel(x, trainSet) ;
		plot(rmodel)  # ,compress=T,branch=0.5);
		text(rmodel) # , pretty=0);
	}
	else if (x$model == "rf" || x$model == "rfNear"){
		if (rfGraphType == "attrEval") {
			imp<-rfAttrEval(x);
			plotRFStats(imp, plotLine=TRUE, myAxes=all.vars(x$formula)[-1]);
		}
		else if (rfGraphType == "attrEvalCluster") {
			#importance by cluster
			impc<-rfAttrEvalClustering(x, trainSet, clustering);
			plotRFMulti(impc$imp, impc$levels, myAxes=all.vars(x$formula)[-1]);
		}
		else if (rfGraphType == "outliers"){
			out<-rfOutliers(x, trainSet);
			plotRFStats(abs(out), cluster=as.character(trainSet[[all.vars(x$formula)[1]]]));
		}
		else if (rfGraphType == "scaling"){
			dis<-rfProximity(x, outProximity=F);
			#get 4 most important components
			space<-spaceScale(dis, 4);
			#  display 1. in 2. component
			subDim<-c(space$points[,1], space$points[,2]);
			dim(subDim)<-c(length(space$points[,1]),2);
			className <- all.vars(x$formula)[1];
			cluster<-trainSet[as.character(className)];
			plotRFStats(subDim, t(cluster));
		}
		else if (rfGraphType == "prototypes"){
			# 10 most typical cases for each class based on predicted class probability
			best<-classPrototypes(x, trainSet, 10);
			vnorm<-varNormalization(x, trainSet[best$prototypes,]);
			plotRFNorm(vnorm, best$cluster, best$levels, 0.15, myHoriz=TRUE, myAxes=all.vars(x$formula)[-1]);
		}
	}
	else {
		warning("The model provided has no visualization.");
	}
	invisible(x)
}
attrEval <- function(formula, data, estimator, costMatrix = NULL, outputNumericSplits=FALSE, ...)
{
	## find the index of estimator
	isRegression <- FALSE ;
	estDsc <- infoCore(what="attrEval");
	estIndex <- match(estimator, estDsc, nomatch=-1);
	if (estIndex == -1) {
		estDscReg <- infoCore(what="attrEvalReg");
		estIndex <- match(estimator, estDscReg, nomatch=-1);
		if (estIndex == -1) 
			stop("Invalid estimator parameter")
		else 
			isRegression = TRUE ;
	}
	
	# check formula or response index or reponse name
	if (inherits(formula,"formula")) {
		dat <- model.frame(formula, data=data, na.action=na.pass)
		trms <- attr(dat,"terms")
		attributes(trms) <- NULL
		formulaExpanded <- as.formula(trms)
	} 
	else {
		if (is.numeric(formula)) {
			if (formula == round(formula)) {# index of response variable
				classIdx <- formula
				className <- names(data)[classIdx]
			}
			else  stop("The first argument must be a formula or prediction column name or prediction column index.")
		}
		else if (is.character(formula)) { # name of response variable
			classIdx <- match(formula, names(data))
			if (length(classIdx) != 1 || is.na(classIdx)) 
				stop("The first argument must be a formula or prediction column name or prediction column index.")
			className <- names(data[classIdx])
		}
		else stop("The first argument must be a formula or prediction column name or prediction column index.")
		
		dat <- data.frame(data[, classIdx], data[, -classIdx, drop=FALSE])
		names(dat)[1] <- className
		# get formula explicitly to allow storage of all terms and their manipulation
		frml <- paste(className, "~",paste(names(dat)[-1], sep="+",collapse="+"),sep="") 
		formulaExpanded <- as.formula(frml)   
	}
	
	if (!isRegression && !inherits(dat[[1]],"factor")) {
		dat[[1]] <- factor(dat[[1]]);
		cat("Changing dependent variable to factor with levels:",levels(dat[[1]]),"\n");
		warning("Possibly this is an error caused by regression formula and classification attribute estimator or vice versa.")
	}
	if (!isRegression && is.null(costMatrix)) {
		class.lev <- levels(dat[[1]]);
		noClasses <- length(class.lev)
		# create and fill uniform costs matrix
		costMatrix <- 1 - diag(noClasses);
	}
	
	aux <- prepare.Data(dat,formulaExpanded,dependent=TRUE,numericAsOrdered=FALSE,skipNAcolumn=TRUE,skipEqualColumn=FALSE);
	discnumvalues <- aux$discnumvalues;
	discdata <- aux$discdata;
	discmap <- aux$discmap;
	numdata <- aux$numdata;
	nummap <- aux$nummap;
	skipmap<-aux$skipmap
	discAttrNames <- dimnames(discdata)[[2]]
	discValCompressed <- aux$disccharvalues
	discValues <- aux$discValues
	numAttrNames <- dimnames(numdata)[[2]]    
	options <- prepare.Options(...);
	
	#check solaris which cannot handle openMP code and force it to use a single thread
#  versionStr <- paste(version,sep="",collapse="")
#  if (grepl("sun",versionStr,fixed=T) || grepl("solaris",versionStr,fixed=TRUE)) {
#	  options[length(options)+1] <- as.character(1)
#	  names(options)[length(options)] <- "maxThreads"
#  }
	
	checkOptionsValues(options) ;
	optRemain <- checkEstimatorOptions(estimator, options, isRegression) ;
	if (length(optRemain) > 0) warning("Unused options:", paste(names(optRemain), collapse=", "));
	if (isRegression) {
		tmp <- .C("estimateCoreReg",
				noInst= aux$noInst,
				noDiscrete = ncol(discdata),
				noDiscreteValues = as.integer(discnumvalues),
				discreteData = as.integer(discdata), # vector of length noInst*noDiscrete, columnwise
				noNumeric = ncol(numdata),
				numericData = as.double(numdata), # vector of length noInst*noNumeric, columnwise
				discAttrNames = as.character(discAttrNames),
				discValNames = as.character(discValCompressed),
				numAttrNames = as.character(numAttrNames),
				numOptions = length(options),
				optionsName = names(options),
				optionsVal = options,
				selEst = estIndex,
				estDisc = double(ncol(discdata)),
				estNum = double(ncol(numdata)),
				splitPointNum = double(ncol(numdata)),
				NAOK=TRUE,
				PACKAGE="CORElearn"
		);
		# assumes length(estNum) == noNumeric, but estNum[1] for predictor is not used
		if (nummap[1] != 1) stop("no dependent variable in prepared regression data"); 
	}
	else {
		tmp <- .C("estimateCore",
				noInst = aux$noInst,
				noDiscrete = ncol(discdata),
				noDiscreteValues = as.integer(discnumvalues),
				discreteData = as.integer(discdata), # vector of length noInst*noDiscrete, columnwise
				noNumeric = ncol(numdata),
				numericData = as.double(numdata), # vector of length noInst*noNumeric, columnwise
				costs = as.double(costMatrix),
				discAttrNames = as.character(discAttrNames),
				discValNames = as.character(discValCompressed),
				numAttrNames = as.character(numAttrNames),            
				numOptions = length(options),
				optionsName = names(options),
				optionsVal = options,
				selEst = estIndex,
				estDisc = double(ncol(discdata)),
				estNum = double(ncol(numdata)),
				splitPointNum = double(ncol(numdata)),
				NAOK=TRUE,
				PACKAGE="CORElearn"
		);
		# assumes length(estDisc) == noDiscrete, but estDist[1] for class is not used
		if (discmap[1] != 1) stop("no class in prepared data"); # for debugging only
	}
	est <- double(length(discmap) + length(nummap)+length(skipmap));
	est[discmap] <- tmp$estDisc;
	est[nummap] <- tmp$estNum;
	names(est)[discmap] <- discAttrNames
	names(est)[nummap] <- numAttrNames   
	
	if (outputNumericSplits) {
		sp <- double(length(discmap) + length(nummap)+length(skipmap));
		sp[nummap] <- tmp$splitPointNum
		names(sp)[nummap] <- numAttrNames
		return( list(attrEval=est[-c(1,skipmap)], splitPointNum=sp[nummap][-1]))
	}
	else { # output only feature evaluations
		return(est[-c(1,skipmap)])
	}
}
rfAttrEval <- function(model) {
	if (! model$model %in% c("rf","rfNear") ) stop("Only random forest model can evaluate attributes with this function.");
	modelID <- model$modelID
	tmp <- .C("rfAttrEval",
			modelID = as.integer(modelID),
			est = double(model$noDiscrete+model$noNumeric),    
			PACKAGE="CORElearn"
	)
	est <- double(length(model$discmap) + length(model$nummap)+length(model$skipmap));
	est[model$discmap] <- tmp$est[1:length(model$discmap)];
	est[model$nummap] <- tmp$est[(length(model$discmap)+1) : (length(model$discmap) + length(model$nummap)) ];
	names(est)[model$discmap] <- model$discAttrNames
	names(est)[model$nummap] <- model$numAttrNames   
	est[-c(1,model$skipmap)];
}

rfOOB <- function(model) {
	if (! model$model %in% c("rf","rfNear") ) 
		stop("Only random forest models can output out of bag performance estimators. Current model is of type ", model$model);
	modelID <- model$modelID
	tmp <- .C("rfOOB",
			modelID = as.integer(modelID),
			oobAccuracy = double(1),  
			oobMargin = double(1),
			oobCorrelation = double(1),
			PACKAGE="CORElearn"
	)
	res<-list(accuracy=tmp$oobAccuracy, margin=tmp$oobMargin, correlation=tmp$oobCorrelation)
	return(res)
}


ordEval <- function(formula, data, file=NULL, rndFile=NULL, variant=c("allNear","attrDist1","classDist1"), ...)
{
	# check formula or response index or reponse name
	if (inherits(formula,"formula")) {
		dat <- model.frame(formula, data=data, na.action=na.pass)
		trms <- attr(dat,"terms")
		attributes(trms) <- NULL
		formulaExpanded <- as.formula(trms)
	} 
	else {
		if (is.numeric(formula)) {
			if (formula == round(formula)) {# index of response variable
				classIdx <- formula
				className <- names(data)[classIdx]
			}
			else  stop("The first argument must be a formula or prediction column name or prediction column index.")
		}
		else if (is.character(formula)) { # name of response variable
			classIdx <- match(formula, names(data))
			if (length(classIdx) != 1 || is.na(classIdx)) 
				stop("The first argument must be a formula or prediction column name or prediction column index.")
			className <- names(data[classIdx])
		}
		else stop("The first argument must be a formula or prediction column name or prediction column index.")
		
		dat <- data.frame(data[, classIdx], data[, -classIdx, drop=FALSE])
		names(dat)[1] <- className
		# get formula explicitly to allow storage of all terms and their manipulation
		frml <- paste(className, "~",paste(names(dat)[-1], sep="+",collapse="+"),sep="") 
		formulaExpanded <- as.formula(frml)   
	}	
	variant <- match.arg(variant)
	variantIdx=match(variant,eval(formals()$variant),nomatch=-1)
	if (!inherits(dat[[1]],"factor")) {
		dat[[1]] <- factor(dat[[1]]);
	}
	class.lev <- levels(dat[[1]]);
	aux <- prepare.Data(dat, formulaExpanded,dependent=TRUE, numericAsOrdered=TRUE,skipNAcolumn=TRUE,skipEqualColumn=TRUE);
	discnumvalues <- aux$discnumvalues;
	discdata <- aux$discdata;
	discmap <- aux$discmap;
	discAttrNames <- dimnames(discdata)[[2]]
	discValNames <- aux$disccharvalues
	options <- prepare.Options(...);
	checkOptionsValues(options) ;
	optRemain <- checkOrdEvalOptions(options)   
	if (length(optRemain) > 0) warning("Unused options:", paste(names(optRemain), collapse=", "));
	noAttr <- ncol(discdata) - 1
	maxAttrValues <- max(discnumvalues[-1])+1	
	statNames<-getStatNames() ;
	noStats <- length(statNames)  ## we get 8 statistics about random normalizers
	tmp <- .C("ordEvalCore",
			noInst =  aux$noInst,
			noDiscrete = ncol(discdata),
			noDiscreteValues = as.integer(discnumvalues),
			discreteData = as.integer(discdata), # vector of length noInst*noDiscrete, columnwise
			discAttrNames = as.character(discAttrNames),
			discValNames = as.character(discValNames),
			numOptions = length(options),
			optionsName = names(options),
			optionsVal = options,
			reinfPos = double(noAttr * maxAttrValues),
			reinfNeg = double(noAttr * maxAttrValues),
			anchor = double(noAttr * maxAttrValues),
			rndReinfPos = double(noAttr * maxAttrValues * noStats),
			rndReinfNeg = double(noAttr * maxAttrValues * noStats),
			rndAnchor = double(noAttr * maxAttrValues * noStats),
			noAV = integer(noAttr * maxAttrValues),
			file = as.character(file),
			rndFile = as.character(rndFile),
			variant = as.integer(variantIdx),
			NAOK=TRUE,
			PACKAGE="CORElearn"
	)
	attrNames <- names(dat)[-1]
	attrMap <- (discmap[-1]) - 1
	attrMapLen <- length(attrMap)
	avNames <- c(1:(maxAttrValues-1),"all")
	avMap <- 1:(maxAttrValues-1)
	avMapLen <- length(avMap)
	reinfPos <- matrix(tmp$reinfPos, nrow=noAttr, ncol=maxAttrValues,dimnames=list(attrNames,avNames));
	reinfNeg <- matrix(tmp$reinfNeg, nrow=noAttr, ncol=maxAttrValues,dimnames=list(attrNames,avNames));
	anchor <- matrix(tmp$anchor, nrow=noAttr, ncol=maxAttrValues,dimnames=list(attrNames,avNames));
	noAV <- matrix(tmp$noAV, nrow=noAttr, ncol=maxAttrValues,dimnames=list(attrNames,avNames));
	rndReinfPos <- array(tmp$rndReinfPos, dim=c(noAttr, maxAttrValues,noStats), dimnames=list(attrNames, avNames, statNames));
	rndReinfNeg <- array(tmp$rndReinfNeg, dim=c(noAttr, maxAttrValues,noStats), dimnames=list(attrNames, avNames, statNames)) ;
	rndAnchor <- array(tmp$rndAnchor, dim=c(noAttr, maxAttrValues,noStats), dimnames=list(attrNames, avNames, statNames));
	rndReinfPosAttr=matrix(rndReinfPos[attrMap,maxAttrValues,],nrow=noAttr, ncol=noStats,dimnames=list(attrNames,statNames))
	rndReinfNegAttr=matrix(rndReinfNeg[attrMap,maxAttrValues,],nrow=noAttr, ncol=noStats,dimnames=list(attrNames,statNames))
	rndAnchorAttr=matrix(rndAnchor[attrMap,maxAttrValues,],nrow=noAttr, ncol=noStats,dimnames=list(attrNames,statNames))
	res<-list(reinfPosAV=reinfPos[attrMap,avMap, drop=FALSE], 
			reinfNegAV=reinfNeg[attrMap,avMap, drop=FALSE], 
			anchorAV=anchor[attrMap,avMap, drop=FALSE], 
			noAV = noAV[attrMap,avMap, drop=FALSE],
			reinfPosAttr=reinfPos[attrMap,maxAttrValues, drop=FALSE], 
			reinfNegAttr=reinfNeg[attrMap,maxAttrValues, drop=FALSE], 
			anchorAttr=anchor[attrMap,maxAttrValues, drop=FALSE],
			noAVattr = noAV[attrMap,maxAttrValues, drop=FALSE], 
			rndReinfPosAV=rndReinfPos[attrMap,avMap, , drop=FALSE], 
			rndReinfNegAV=rndReinfNeg[attrMap,avMap, , drop=FALSE], 
			rndAnchorAV=rndAnchor[attrMap,avMap, , drop=FALSE],
			rndReinfPosAttr=rndReinfPosAttr,
			rndReinfNegAttr=rndReinfNegAttr, 
			rndAnchorAttr=rndAnchorAttr,
			attrNames= attrNames, 
			valueNames=aux$discValues[-1], 
			noAttr=length(attrNames),
			ordVal=maxAttrValues-1,
			variant=variant,
			file=file, 
			rndFile=rndFile,
			formula=aux$formulaOut
	);
	class(res) <- "ordEval"  
	return(res)
}
plotInstEval<-function(oeInstFile, oeInstRndFile,  noAttr, ...) {
	inst<-read.table(oeInstFile,header=FALSE,sep=",",colClasses="character",strip.white=TRUE,na.strings=c("NA","?"))
	instNorm<-read.table(oeInstRndFile,header=FALSE,sep=",",colClasses="character",strip.white=TRUE,na.strings=c("NA","?"))
	#noAttr <- length(ordEvalData$noAVattr)
	#ordVal <- ncol(ordEvalData$reinfPosAV)
	statNames<-getStatNames() ;
	noStats <- length(statNames)
	ord<-list()
	noInst <- nrow(inst)/(noAttr+1)
	for (i in 1:noInst) {
		className <- as.character(trimSpaces(inst[(i-1)*(noAttr+1)+1, 1]))
		classValue <- as.character(trimSpaces(inst[(i-1)*(noAttr+1)+1, 2]))
		attrName <- as.character(trimSpaces(inst[((i-1)*(noAttr+1)+2):(i*(noAttr+1)), 1]))
		valueName <- as.character(trimSpaces(inst[((i-1)*(noAttr+1)+2):(i*(noAttr+1)), 2]))
		reinfPos <- as.numeric(inst[((i-1)*(noAttr+1)+2):(i*(noAttr+1)), 3])
		reinfNeg <- as.numeric(inst[((i-1)*(noAttr+1)+2):(i*(noAttr+1)), 4])
		anchor <- as.numeric(inst[((i-1)*(noAttr+1)+2):(i*(noAttr+1)), 5])
		rndReinfPos<-list()
		rndReinfNeg<-list()
		rndAnchor<-list()
		for (iA in 1:noAttr){
			rndReinfPos[[iA]]<-as.numeric(instNorm[(i-1)*(noAttr+1)+1+iA,2:(1+noStats)])
			rndReinfNeg[[iA]]<-as.numeric(instNorm[(i-1)*(noAttr+1)+1+iA,(2+noStats):(1+2*noStats)])
			rndAnchor[[iA]]<-as.numeric(instNorm[(i-1)*(noAttr+1)+1+iA,(2+2*noStats):(1+3*noStats)])
		}
		ord[[i]]<-list(className=className,classValue=classValue,attributeName=attrName,valueName=valueName,reinfPos=reinfPos,reinfNeg=reinfNeg,anchor=anchor,
				rndReinfPos=rndReinfPos,rndReinfNeg=rndReinfNeg,rndAnchor=rndAnchor)
	} 
	oeInst(ord, noAttr, ...)
}
plotOrdEval<-function(file, rndFile=NULL, ...){
	# read data from files and transform the two tables to internal object as returned by ordEval
	ord<-read.table(file,header=TRUE,sep=",",strip.white=TRUE)
	if (!is.null(rndFile))
		ordNorm<-read.table(rndFile,header=TRUE,sep=",",strip.white=TRUE)
	## extract number of attributes and values from first column
	name <- ord[,1]
	dup <- duplicated(name)
	for (i in 2:length(dup))
		if (dup[i])
			break ;
	ordVal <- i-3
	noAttr <- length(unique(name)) - ordVal
	statNames<-getStatNames()
	noStats <- length(statNames)  ## we get 8 statistics about random normalizers
	attrNames <- c()
	avNames <- c(1:ordVal)
	for (iA in 1:noAttr) {
		attrNames[iA] <- as.character(ord[(iA-1)*(ordVal+1)+1,1]) 
	}
	reinfPosAV <- matrix(0, nrow=noAttr, ncol=ordVal,dimnames=list(attrNames,avNames));
	reinfNegAV <- matrix(0, nrow=noAttr, ncol=ordVal,dimnames=list(attrNames,avNames));
	anchorAV <- matrix(0, nrow=noAttr, ncol=ordVal,dimnames=list(attrNames,avNames));
	noAV <- matrix(0, nrow=noAttr, ncol=ordVal,dimnames=list(attrNames,avNames));
	reinfPosAttr <- array(0,dim=c(noAttr),dimnames=list(attrNames));
	reinfNegAttr <- array(0, dim=c(noAttr),dimnames=list(attrNames));
	anchorAttr <- array(0, dim=c(noAttr),dimnames=list(attrNames));
	noAVattr <- array(0, dim=c(noAttr),dimnames=list(attrNames));
	rndReinfPosAV <- array(0, dim=c(noAttr, ordVal,noStats), dimnames=list(attrNames, avNames, statNames));
	rndReinfNegAV <- array(0, dim=c(noAttr, ordVal,noStats), dimnames=list(attrNames, avNames, statNames)) ;
	rndAnchorAV <- array(0, dim=c(noAttr, ordVal,noStats), dimnames=list(attrNames, avNames, statNames));
	rndReinfPosAttr <- array(0, dim=c(noAttr, noStats), dimnames=list(attrNames,  statNames));
	rndReinfNegAttr <- array(0, dim=c(noAttr, noStats), dimnames=list(attrNames,  statNames)) ;
	rndAnchorAttr <- array(0, dim=c(noAttr,noStats), dimnames=list(attrNames,  statNames));
	valueNames<-list()
	for (iA in 1:noAttr) {
		#attrNames[iA] <- ord[(iA-1)*(ordVal+1)+1,1]
		valueNames[[iA]] <- as.character(ord[(2+(iA-1)*(ordVal+1)):(iA*(ordVal+1)),1])
		noAV[iA,] <- ord[(2+(iA-1)*(ordVal+1)):(iA*(ordVal+1)),5]       
		for(i in 1:ordVal) {
			reinfPosAV[iA,i]  <- ord[(iA - 1) * (ordVal + 1) + i + 1, 2]
			reinfNegAV[iA,i] <- ord[(iA - 1) * (ordVal + 1) + i + 1, 3]
			anchorAV[iA,i]   <- ord[(iA - 1) * (ordVal + 1) + i + 1, 4]
			
			if (!is.null(rndFile)) {
				rndReinfPosAV[iA,i,] <- as.numeric(ordNorm[(iA - 1) * (ordVal + 1) + i + 1, 2:(noStats+1)])
				rndReinfNegAV[iA,i,]  <- as.numeric(ordNorm[(iA - 1) * (ordVal + 1) + i + 1, (2+noStats):(1+2*noStats)])
				rndAnchorAV[iA,i,]      <- as.numeric(ordNorm[(iA - 1) * (ordVal + 1) + i + 1, (2+2*noStats):(1+3*noStats)])
			}
		}
		i <- 0
		
		noAVattr[iA] <- ord[(1+(iA-1)*(ordVal+1)),5]               
		
		reinfPosAttr[iA]  <- ord[(iA - 1) * (ordVal + 1) + i + 1, 2]
		reinfNegAttr[iA] <- ord[(iA - 1) * (ordVal + 1) + i + 1, 3]
		anchorAttr[iA]   <- ord[(iA - 1) * (ordVal + 1) + i + 1, 4]
		
		if (!is.null(rndFile)) {
			rndReinfPosAttr[iA,] <- as.numeric(ordNorm[(iA - 1) * (ordVal + 1) + i + 1, 2:(noStats+1)])
			rndReinfNegAttr[iA,]  <- as.numeric(ordNorm[(iA - 1) * (ordVal + 1) + i + 1, (2+noStats):(1+2*noStats)])
			rndAnchorAttr[iA,]      <- as.numeric(ordNorm[(iA - 1) * (ordVal + 1) + i + 1, (2+2*noStats):(1+3*noStats)])
		}
	}
	oeObj <- list(reinfPosAV=reinfPosAV, reinfNegAV=reinfNegAV, anchorAV=anchorAV, noAV = noAV,
			reinfPosAttr=reinfPosAttr, reinfNegAttr=reinfNegAttr, anchorAttr=anchorAttr, noAVattr = noAVattr,
			rndReinfPosAV=rndReinfPosAV, rndReinfNegAV=rndReinfNegAV, rndAnchorAV=rndAnchorAV,
			rndReinfPosAttr=rndReinfPosAttr, rndReinfNegAttr=rndReinfNegAttr, rndAnchorAttr=rndAnchorAttr,
			attrNames= attrNames, valueNames=valueNames, noAttr=length(attrNames),ordVal=ordVal,variant=NULL,file=file, rndFile=rndFile      
	) 
	class(oeObj) <- "ordEval"
	plot(oeObj,...)  ## call of plot.ordEval
}

plot.ordEval<-function(x, graphType=c("avBar", "attrBar", "avSlope"), ...) {
	
	graphType<-match.arg(graphType)
	
	if (graphType=="avSlope")
		avSlopeObject(x,  ...)
	else if (graphType=="avBar" )
		avNormBarObject(x, ...)
	else if (graphType=="attrBar")
		attrNormBarObject(x, ...)
	invisible(x)
}

printOrdEval<-function(x) {
	object <- x
	maxAttrChars <- max(nchar(c(object$attrNames,"Attribute")))
	maxAVChars <- max(nchar(c(unlist(object$valueNames),"Value")))
	header <- paste(sprintf("%*s %*s",maxAttrChars,"Attribute",maxAVChars,"Value"),
			sprintf("%6s %6s %6s %6s  %6s %6s %6s %6s  %6s %6s %6s %6s","Down","Down_p","Down_l","Down_h","Up","Up_p","Up_l","Up_h","Anchor","Anch_p","Anch_l","Anch_h")
			,sep=" ")
	cat(header,"\n")
	for (a in 1:object$noAttr) {
		line <- paste(sprintf("%*s %*s",maxAttrChars,object$attrNames[a],maxAVChars,"all"),
				sprintf("%6.4f %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f %6.4f", 
						object$reinfPosAttr[a], 
						object$rndReinfPosAttr[a,"p-value"],
						object$rndReinfPosAttr[a,"lowPercentile"],
						object$rndReinfPosAttr[a,"highPercentile"],
						object$reinfNegAttr[a],
						object$rndReinfNegAttr[a,"p-value"],
						object$rndReinfNegAttr[a,"lowPercentile"],
						object$rndReinfNegAttr[a,"highPercentile"],
						object$anchorAttr[a],
						object$rndAnchorAttr[a,"p-value"],
						object$rndAnchorAttr[a,"lowPercentile"],
						object$rndAnchorAttr[a,"highPercentile"]),               
				sep="  ") ;
		cat(line,"\n")
		
		for (v in 1:object$ordVal){
			line <- paste(sprintf("%*s %*s",maxAttrChars," ",maxAVChars,object$valueNames[[a]][v]),
					sprintf("%6.4f %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f %6.4f", 
							object$reinfPosAV[a,v],
							object$rndReinfPosAV[a,v,"p-value"],
							object$rndReinfPosAV[a,v,"lowPercentile"],
							object$rndReinfPosAV[a,v,"highPercentile"],
							object$reinfNegAV[a,v],
							object$rndReinfNegAV[a,v,"p-value"],
							object$rndReinfNegAV[a,v,"lowPercentile"],
							object$rndReinfNegAV[a,v,"highPercentile"],
							object$anchorAV[a,v],
							object$rndAnchorAV[a,v,"p-value"],
							object$rndAnchorAV[a,v,"lowPercentile"],
							object$rndAnchorAV[a,v,"highPercentile"]),            
					sep="  ")
			cat(line,"\n")
		}
	}
}


modelEval <- function(model=NULL, correctClass, predictedClass, predictedProb=NULL, costMatrix=NULL, priorClProb = NULL, avgTrainPrediction = NULL, beta=1) {
	if (is.null(predictedClass) && is.null(predictedProb)) {
		warning("Only one of the predictedClass and predictedProb parameters can be NULL")
		return(NULL) ;
	}      
	if (is.null(model)) {
		if (!is.null(avgTrainPrediction)) {
			return(modelEvaluationReg.Core(correctClass,predictedClass,avgTrainPrediction))
		}
		else { 
			if (is.null(predictedClass))
				predictedClass <- levels(correctClass)[apply(predictedProb, 1, which.max)]            
			return(modelEvaluationClass.Core(correctClass,predictedClass,predictedProb,costMatrix,priorClProb,beta))
		}
	}
	if (class(model) != "CoreModel"){
		warning("Only models of type CoreModel can be evaluated with this type of call. Others shall supply NULL for parameter model, and provide value of avgTrainPrediction in case of regression.")
		return(NULL) ;
	}
	if (model$model == "regTree") {
		if (is.null(avgTrainPrediction))
			avgTrainPrediction <- model$avgTrainPrediction
		return(modelEvaluationReg.Core(correctClass,predictedClass,model$avgTrainPrediction))
	}
	else {
		if (is.null(priorClProb))
			priorClProb <- model$priorClassProb
		if (is.null(predictedClass))
			predictedClass <- model$class.lev[apply(predictedProb, 1, which.max)]
		return(modelEvaluationClass.Core(correctClass,predictedClass,predictedProb,costMatrix,priorClProb,beta))
	}
}  

modelEvaluationClass.Core <- function(correctClass, predictedClass, predictedProb=NULL, costMatrix=NULL, priorClassProb=NULL, beta=1) {
	# common vector of levels
	if (inherits(correctClass,"factor")) {
		levelsCorrect<-levels(correctClass)
	} else {
		levelsCorrect<-sort(unique(correctClass))
	}
	if (inherits(predictedClass,"factor")) {
		levelsPredicted<-levels(predictedClass)
	} else {
		levelsPredicted<-sort(unique(predictedClass))
	}
	levelsBoth<-union(levelsCorrect,levelsPredicted)
	# some data validity checks
	correctClass<-factor(correctClass,levels=levelsBoth)
	if (any(is.na(correctClass)))
		stop("Correct class should not contain NA values.")
	noClasses <- length(levelsBoth)
	predictedClass<-factor(predictedClass,levels=levelsBoth)
	if (any(is.na(predictedClass)))
		stop("Predicted class should not contain NA values.")
	noInst <- length(correctClass)
	if (is.null(predictedProb)){
		## create and fill the prediction matrix
		predictedProb <- matrix(0, nrow=noInst, ncol=noClasses)
		for (i in 1:noInst)
			predictedProb[i, predictedClass[i]] <- 1
	}
	if (is.null(costMatrix)) {
		## create and fill uniform costs matrix
		costMatrix <- 1 - diag(noClasses)
	}
	if (is.null(priorClassProb))
		priorClassProb <- table(correctClass)/noInst
	tmp <- .C("modelEvaluate",
			noInst = length(correctClass),
			correctClass = as.integer(correctClass),
			# predictedClass = as.integer(predictedClass), # computed from predictedProb and CostMatrix
			predictedProb = as.double(predictedProb),
			costMatrix = as.double(costMatrix),
			noClasses = as.integer(noClasses), 
			priorClassProb = as.double(priorClassProb),
			accuracy = double(1),
			avgCost = double(1),
			infScore = double(1),
			auc = double(1),
			predMatrix = integer(noClasses * noClasses),
			sensitivity = double(1),
			specificity = double(1),
			brier = double(1),
			kappa = double(1), 
			precision = double(1),
			Gmean = double(1),    
			KS = double(1),
			TPR = double(1),
			FPR = double(1),
			NAOK=TRUE,
			PACKAGE="CORElearn"
	)
	recall = tmp$sensitivity
	denominator <- (beta*beta * recall + tmp$precision)
	if (denominator == 0)
		Fmeasure <- 0
	else
		Fmeasure = (1+beta*beta)*recall*tmp$precision / denominator
	predMx <- matrix(tmp$predMatrix, nrow = noClasses, ncol=noClasses, dimnames = list(levels(correctClass),levels(correctClass)))
	list(accuracy = tmp$accuracy, averageCost = tmp$avgCost, informationScore = tmp$infScore,
			AUC = tmp$auc, predictionMatrix = predMx, sensitivity = tmp$sensitivity,
			specificity = tmp$specificity, brierScore = tmp$brier, kappa = tmp$kappa,
			precision = tmp$precision, recall = tmp$sensitivity, Fmeasure = Fmeasure, 
			Gmean = tmp$Gmean, KS = tmp$KS, TPR = tmp$TPR, FPR = tmp$FPR)
}

modelEvaluationReg.Core <- function(correct, predicted, avgTrainPredicted) {
	noInst <- length(correct) ;
	tmp <- .C("modelEvaluateReg",
			noInst = length(correct),
			correct = as.double(correct),
			predicted = as.double(predicted),
			avgPredicted = as.double(avgTrainPredicted), 
			MSE = double(1),
			RMSE = double(1),
			MAE = double(1),
			RMAE = double(1),
			NAOK=TRUE,
			PACKAGE="CORElearn"
	)
	list(MSE = tmp$MSE, RMSE = tmp$RMSE, MAE = tmp$MAE, RMAE = tmp$RMAE)
}

paramCoreIO <- function(model, fileName, io=c("read","write")) {
	io = match.arg(io)
	tmp <- .C("optionsInOut",
			modelID = as.integer(model$modelID),
			fileName=as.character(fileName),
			io=as.character(io),
			NAOK=FALSE,
			PACKAGE="CORElearn"
	)
	invisible(tmp)
}

saveRF <- function(model, fileName) {
	if (model$model != "rf") stop("Only random forest model can be saved at the moment.");
	modelID <- model$modelID
	tmp <- .C("saveRF",
			modelID = as.integer(modelID),
			fileName=as.character(fileName),
			PACKAGE="CORElearn"
	)
	invisible(tmp)
}

loadRF <- function(formula, data, fileName) {
	model="rf"
	# check formula or response index or reponse name
	if (inherits(formula,"formula")) {
		dat <- model.frame(formula, data=data, na.action=na.pass)
		trms <- attr(dat,"terms")
		attributes(trms) <- NULL
		formulaExpanded <- as.formula(trms)
	} 
	else {
		if (is.numeric(formula)) {
			if (formula == round(formula)) {# index of response variable
				classIdx <- formula
				className <- names(data)[classIdx]
			}
			else  stop("The first argument must be a formula or prediction column name or prediction column index.")
		}
		else if (is.character(formula)) { # name of response variable
			classIdx <- match(formula, names(data))
			if (length(classIdx) != 1 || is.na(classIdx)) 
				stop("The first argument must be a formula or prediction column name or prediction column index.")
			className <- names(data[classIdx])
		}
		else stop("The first argument must be a formula or prediction column name or prediction column index.")
		
		dat <- data.frame(data[, classIdx], data[, -classIdx, drop=FALSE])
		names(dat)[1] <- className
		# get formula explicitly to allow storage of all terms and their manipulation
		frml <- paste(className, "~",paste(names(dat)[-1], sep="+",collapse="+"),sep="") 
		formulaExpanded <- as.formula(frml)   
	}
	if (!inherits(dat[[1]],"factor")) {
		dat[[1]] <- factor(dat[[1]]);
		cat("Changing dependent variable to factor with levels:",levels(dat[[1]]),"\n");
	}
	class.lev <- levels(dat[[1]]);
	noClasses <- length(class.lev);
	
	tmp <- .C("readRF",
			fileName=as.character(fileName),
			modelID = integer(1),
			PACKAGE="CORElearn"
	);
	if (tmp$modelID == -1) {
		return(NULL)
	}
	res <- list(modelID=tmp$modelID, class.lev=class.lev, model=model, formula=formula,  noClasses = noClasses)
	class(res) <- "CoreModel"
	res
}

getRFsizes <- function(model, type=c("size", "sumdepth")) {
	if (model$model != "rf") stop("The model must be a random forest.");
	type <- match.arg(type)
	switch(type,
			size=.Call("exportSizesRF", as.integer(model$modelID), PACKAGE="CORElearn"),
			sumdepth=.Call("exportSumOverLeavesRF", as.integer(model$modelID), PACKAGE="CORElearn"))
}   

getCoreModel <- function(model) {
	if (model$model != "rf") stop("The model must be a random forest.");
	.Call("exportModel", as.integer(model$modelID), PACKAGE="CORElearn")    
}

calibrate <- function(correctClass, predictedProb, class1=1, method = c("isoReg","binIsoReg","binning","mdlMerge"), weight=NULL,noBins=10, assumeProbabilities=FALSE){
	noClasses <- length(levels(correctClass)) ;
	method <- match.arg(method)
	methodIdx = match(method, eval(formals()$method), nomatch=-1)
	if (assumeProbabilities==TRUE && any(predictedProb >1.0 | predictedProb<0))
		stop("Predicted probabilities in predictedValues are expected to be in [0,1] range.")
	noInst <- length(correctClass) ;
	if (is.null(weight)) {
		weight<-numeric(noInst)
		weight[]<-1
	}
	# class1 can be either class name (factor) or its index
	if (is.factor(class1))
		class1idx<-match(class1, levels(correctClass),nomatch=-1)
	else {
		class1idx<-class1
		class1<-factor(levels(correctClass)[class1idx],levels=levels(correctClass))       
	}
	# convert true class to a vector of 0 and 1
	tc<-integer(length(correctClass))
	tc[]<-0
	tc[correctClass==class1]<-1
	
	tmp <- .C("calibrate",
			methodIdx = as.integer(methodIdx),
			noInst = as.integer(noInst),
			correctClass = as.integer(tc),
			predictedProb = as.double(predictedProb),
			weight=as.double(weight),
			noBins=as.integer(noBins),
			noIntervals = integer(1),
			interval = double(noInst),
			calProb = double(noInst),
			NAOK=TRUE,
			PACKAGE="CORElearn"
	)
	if (assumeProbabilities == TRUE)
		tmp$interval[tmp$noIntervals] <- 1 # set sentinel for probabilities
	list(interval = tmp$interval[1:tmp$noIntervals], calProb = tmp$calProb[1:tmp$noIntervals])
}

applyDiscretization <- function(data, boundsList, noDecimalsInValueName=2) {
	if (is.null(boundsList))
		return(data)
	
	for (i in 1:length(boundsList)) {
		attrName <- names(boundsList)[i]
		discValues <- apply( outer(data[,attrName], boundsList[[i]], ">"), 1, sum) 
		data[,attrName] <- factor(discValues, levels=0:length(boundsList[[i]]))
		repeat {
			levels(data[,attrName]) <- intervalNames(boundsList[[i]], noDecimalsInValueName)
			if (length(unique(levels(data[,attrName])))==length(boundsList[[i]])+1)
				break
			else
				noDecimalsInValueName <- noDecimalsInValueName +1	
		}
	}
	data
}


discretize <- function(formula, data, method=c("greedy", "equalFrequency", "equalWidth"), estimator=NULL, 
		discretizationLookahead=3,discretizationSample=0, maxBins=0, equalDiscBins=4, ...)

{
	method <- match.arg(method)
	methodIdx = match(method,as.vector(formals()$"method",mode="character")[-1],nomatch=-1)
	
	isRegression <- NULL
	
	# check formula or response index or reponse name
	if (inherits(formula,"formula")) {
		dat <- model.frame(formula, data=data, na.action=na.pass)
		trms <- attr(dat,"terms")
		attributes(trms) <- NULL
		formulaExpanded <- as.formula(trms)
	} 
	else {
		if (is.numeric(formula)) {
			if (formula == round(formula)) {# index of response variable
				classIdx <- formula
				className <- names(data)[classIdx]
			}
			else  stop("The first argument must be a formula or prediction column name or prediction column index.")
		}
		else if (is.character(formula)) { # name of response variable
			classIdx <- match(formula, names(data))
			if (length(classIdx) != 1 || is.na(classIdx)) 
				stop("The first argument must be a formula or prediction column name or prediction column index.")
			className <- names(data[classIdx])
		}
		else stop("The first argument must be a formula or prediction column name or prediction column index.")
		
		dat <- data.frame(data[, classIdx], data[, -classIdx, drop=FALSE])
		names(dat)[1] <- className
		# get formula explicitly to allow storage of all terms and their manipulation
		frml <- paste(className, "~",paste(names(dat)[-1], sep="+",collapse="+"),sep="") 
		formulaExpanded <- as.formula(frml)   
	}
	
	## find the index of estimator
	if (is.null(estimator)) {
		isRegression <- ! inherits(dat[[1]], "factor")
		if (isRegression) # regression
			estimator <- "RReliefFexpRank"
		else
			estimator <- "ReliefFexpRank"		
	} 		
	estDsc <- infoCore(what="attrEval");
	estIndex <- match(estimator, estDsc, nomatch=-1);
	if (estIndex == -1) {
		estDscReg <- infoCore(what="attrEvalReg");
		estIndex <- match(estimator, estDscReg, nomatch=-1);
		if (estIndex == -1) 
			stop("Invalid estimator parameter")
		else isRegression <- TRUE
	}
	else isRegression <- FALSE
	
	if (method == "greedy") {
		if (!is.numeric(maxBins))
			stop("The maximal number of bins shall be an integer or integer vector of length equal to the number of numeric attributes.")
		if (any(maxBins<0 | maxBins==1))
			stop("The maximal number of bins shall be 0 (don't care) or an integer >=2")	
	}
	else { # if (method == "equalFrequency" || method="equalWidth" ) 
		if (!is.numeric(equalDiscBins))
			stop("The number of bins (equalDiscBins) shall be an integer or an integer vector of length equal to the number of numeric attributes.")
		if (any(equalDiscBins<2))
			stop("The number of bins (equalDiscBins) shall be an integer >=2")
	}
	
	
	#if (is.null(isRegression)) # in case of equal width or equal frequency discretization 
	#	isRegression <- ! inherits(dat[[1]], "factor")
	
	if (!isRegression && !inherits(dat[[1]],"factor")) {
		dat[[1]] <- factor(dat[[1]]);
		cat("Changing dependent variable to factor with levels:", levels(dat[[1]]),"\n");
		warning("Possibly this is an error caused by regression formula and classification attribute estimator or vice versa.")
	}
	
	aux <- prepare.Data(dat,formulaExpanded,dependent=TRUE,numericAsOrdered=FALSE,orderedAsNumeric=FALSE, skipNAcolumn=TRUE,skipEqualColumn=FALSE);
	discnumvalues <- aux$discnumvalues;
	discdata <- aux$discdata;
	discmap <- aux$discmap;
	numdata <- aux$numdata;
	nummap <- aux$nummap;
	skipmap<-aux$skipmap
	discAttrNames <- dimnames(discdata)[[2]]
	discValCompressed <- aux$disccharvalues
	discValues <- aux$discValues
	numAttrNames <- dimnames(numdata)[[2]]    
	bounds <- matrix(0, nrow=nrow(numdata),ncol=ncol(numdata))
	options <- prepare.Options(...);
	options[[length(options)+1]] <- discretizationLookahead
	options[[length(options)+1]] <- discretizationSample
	names(options)[(length(options)-1):length(options)] <- c("discretizationLookahead","discretizationSample")
	checkOptionsValues(options) ;
	if (isRegression) {
		if (nummap[1] != 1) stop("No dependent variable in prepared regression data.");
		attr2Disc <- ncol(numdata)-1
	}
	else {
		if (discmap[1] != 1) stop("No class in prepared data."); 
		attr2Disc <- ncol(numdata)
	}
	
	if (attr2Disc==0)
		return(NULL)

	if (method == "greedy") 
		maxBins <- rep(maxBins, length.out=attr2Disc)
	else 
	  maxBins <- rep(equalDiscBins, length.out
	                 =attr2Disc)	
	
	tmp <- .C("discretize",
			methodIdx = as.integer(methodIdx),
			isRegression = as.integer(isRegression),
			noInst = aux$noInst,
			noDiscrete = ncol(discdata),
			noDiscreteValues = as.integer(discnumvalues),
			discreteData = as.integer(discdata), # vector of length noInst*noDiscrete, columnwise
			noNumeric = ncol(numdata),
			numericData = as.double(numdata), # vector of length noInst*noNumeric, columnwise
			discAttrNames = as.character(discAttrNames),
			discValNames = as.character(discValCompressed),
			numAttrNames = as.character(numAttrNames),            
			numOptions = length(options),
			optionsName = names(options),
			optionsVal = options,
			selEst = estIndex,
			maxBins = as.integer(maxBins),
			noBounds = integer(ncol(bounds)),
			bounds = as.double(bounds), # vector of length noInst*noNumeric, columnwise
			NAOK=TRUE,
			PACKAGE="CORElearn"
	)
	boundsMx <- matrix(tmp$bounds, nrow=nrow(bounds),ncol=ncol(bounds),byrow=FALSE)
	outBounds <- list()
	for (i in 1:ncol(boundsMx)) {
		if (tmp$noBounds[i]>0)
			outBounds[[i]] <-  boundsMx[1:tmp$noBounds[i], i]
		else
			outBounds[[i]] <- NA
	}
	names(outBounds) <- numAttrNames
	
	if (isRegression)
		outBounds[[1]] <- NULL
	return(outBounds)
}

noEqualRows <- function(data1, data2, tolerance=1e-5, countOnce=TRUE) {
	if (ncol(data1) != ncol(data2))
		stop("Only data sets with equal number of columns can be compared.")
	d1 <- data.matrix(data1)
	d2 <- data.matrix(data2)
	replaceNA <- (max(d1,d2, na.rm=TRUE)+tolerance)*1.0001 # larger value than any existing
	d1[is.na(d1)] <- replaceNA
	d2[is.na(d2)] <- replaceNA
	storage.mode(d1) <- "double"
	storage.mode(d2) <- "double" 
	.Call("noEqualRows", d1, d2, as.integer(nrow(d1)), as.integer(nrow(d2)), as.integer(ncol(d1)), 
			as.double(tolerance), as.integer(countOnce), PACKAGE="CORElearn")
}