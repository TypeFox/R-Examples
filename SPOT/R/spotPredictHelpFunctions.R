###################################################################################
#' Initialize Predictor
#'  
#' Initialization function that should be called at the beginning of any call to the surrogate models (e.g. \code{spotPredictRandomForest}, \code{spotPredictDace}, ...).
#' 
#' @param design matrix, dataframe or vector that contains an experimental design of one or more observations and variables.
#' @param designType specify "data.frame" to force the design into a dataframe. specify "as.matrix" to force a matrix.
#' @param roi region of interest, a data frame or matrix that has a row for each parameter. It describes the region in which the design is chosen. Row names are mandatory, and describe the parameter names.
#' @param pkgs package(s) that need to be loaded. Can be vector of several packages. If this is NULL, nothing is loaded (default).
#' @param name name of the predictor, used only for a startup message. 
#' @param verbosity this specifies whether startup message is produced \code{(verbosity >= 3)} or not \code{(verbosity < 3)}.
#'
#' @return returns the \code{design}, with the chosen format (matrix, dataframe) and transposed if number of columns do not fit to the \code{roi}.
#'
#' @export
#' @keywords internal
###################################################################################
spotInitializePredictor <- function (design,designType,roi,pkgs=NULL,name,verbosity){
	spotWriteLines(verbosity,3,paste(name,"started"));
	if(!is.null(pkgs))
		spotInstAndLoadPackages(pkgs)
	if(is.null(design)) return(NULL)
	if(designType=="data.frame"){
		design<-data.frame(design)
	}else if(designType=="as.matrix"){
		design<-as.matrix(design)
	}
	if(ncol(design)!=nrow(roi)){ #ugly, but necessary because R is strange about converting one dimensional vectors to dataframes (confusing rows/columns)
		design<-t(design)
		#bugfix, because t() returns a matrix.... >_<
		if(designType=="data.frame"){	
			design<-data.frame(design)
		}		
	}
	colnames(design)<-row.names(roi)
	rownames(design)<-NULL
	design
}

###################################################################################
#' Force into Factors: X
#'  
#' Forces the type of parameters in the sample locations x to factors, if type[i] is specified to "FACTOR"
#' 
#' @param x sample locations
#' @param type string of type information, factors should be  "FACTOR", rest is ignored
#' @param xNames names of the parameters in x
#'
#' @return returns a list with component \code{x} (data frame x, now with factors) \code{levels} (list of levels of the corresponding factors).
#'
#' @keywords internal
###################################################################################
spotForceFactorX<- function (x,type,xNames){
	lvls <- list()
	for(i in 1:ncol(x)){ 
		if(type[i]=="FACTOR"){
			if(!is.factor(x[[xNames[i]]]))
				x[[xNames[i]]]<-as.factor(round(x[[xNames[i]]]))
			lvls[[i]] <- levels(x[[xNames[i]]])
		}	
	}
	list(x=x,levels=lvls)
}


###################################################################################
#' Force into Factors: Design
#'  
#' Forces the type of parameters in the design to be predicted to factors, if type[i] is specified to "FACTOR"
#' 
#' @param design prediction locations
#' @param type string of type information, factors should be  "FACTOR", rest is ignored
#' @param xNames names of the parameters in x
#' @param lvls levels to be used when forcing the parameter values to factor type
#'
#' @return returns the design variable, but now with factors where necessary.
#'
#' @keywords internal
###################################################################################
spotForceFactorDesign<- function (design,type,xNames,lvls){
	for(i in 1:ncol(design)){ #designs produce numerics, but earth can deal with factors, thus force factor columns "as.factor"
		if(type[i]=="FACTOR"){
			if(!is.factor(design[[xNames[i]]]))
				design[[xNames[i]]]<-factor(round(design[[xNames[i]]]),levels=lvls[[i]])
		}	
	}
	design
}

###################################################################################
#' Create Formula for Coding
#'  
#' Create the formulas used by spotPredictLm and similar functions to code input data.
#' Uncoded data has column names as chosen in ROI, coded data iterates names  "x1", "x2",... .
#' 
#' @param X sample locations to be coded
#' @param name names of uncoded data columns
#' @param type data type of each data column. Those chosen with "FACTOR" will not be coded. 
#'
#' @return returns a list of formulas. One formula for each column in X which is not a "FACTOR".
#'
#' @keywords internal
###################################################################################
spotCodeFormulas <- function (X,name,roi){
	type <- roi$type
	fmla = NULL
	for (i in 1:length(type)){
		if(type[i]!="FACTOR"){ #only non-factors should be coded. factors can not be coded, as they have no order.
			a <- roi$lower[i]
			b <- roi$upper[i]
			v1 <- mean(c(a,b))
			v2 <- (b-a)/2	
			fmla <- c(fmla,as.formula(paste(paste("x",i,sep=""),"~ (", name[i], " - ", v1, ")/", v2)))		
		}	
	}
	fmla
}

###################################################################################
#' Model Training Interface
#'  
#' This function is used to interface the \code{spotPredict*} (e.g.: \code{\link{spotPredictRandomForest}}) functions.
#' It serves as a simple interface to train surrogate models, as need by the user.
#' 
#' @param X Matrix of sample locations. Rows for points, columns for variables.
#' @param y Observations at sample locations (vector)
#' @param method which model function should be used specified as a string, e.g., "spotPredictRandomForest", "spotPredictLm", "spotPredictForrester". Most other "spotPredict*" functions provided by SPOT should work.
#' @param settings additional settings, as used by the corresponding \code{spotPredict*} function. NULL means default model settings are used.
#'
#' @return the list of settings, including the model fit. This list can be employed as the "fit" parameter in \code{\link{spotModel.predict}}.
#'
#' @seealso instead of building a fit to be evaluated by \code{\link{spotModel.predict}} the \code{\link{spotModel.func}} function
#' directly produces functions, which can be evaluated with new data. 
#'
#' @export
#'
#' @examples
#' ## simple test function
#' sphereFunction <- function(X) rowSums(X^2)
#' ## sample locations
#' Xtrain <- matrix(runif(20),,2)
#' ## evaluate on test function
#' y <- sphereFunction(Xtrain) 
#' ## train model
#' model <- spotModel.train(Xtrain,y,"spotPredictLm")
#' ## print model fit
#' print(model$seq.modelFit)
###################################################################################
spotModel.train <- function(X,y,method="spotPredictTree",settings=NULL) {
	X <- data.frame(X)
	xNames<-colnames(X)
	y <- data.frame(y)
	yNames <- paste("y.",colnames(y),sep="")
	colnames(y)<-yNames
	rawB <- cbind(X,y)	
	spotConfig <- list(
		io.verbosity=0
		,seq.predictionModel.func = method
		,alg.roi = spotROI(apply(X,2,min),apply(X,2,max),varnames=xNames)
		,alg.resultColumn = yNames
		,spot.fileMode = FALSE
		,seq.model.variance = FALSE
		)
	spotConfig <- append(settings,spotConfig); 
	spotConfig <- spotConfig[!duplicated(names(spotConfig))]	
	eval(call(method, rawB, NULL, NULL, spotConfig))
}

###################################################################################
#' Model Prediction Interface
#'  
#' This function is used to interface the \code{spotPredict*} (e.g.: \code{\link{spotPredictRandomForest}}) functions.
#' It serves as a simple interface to predict new data with on a surrogate models.
#' 
#' @param X Matrix of locations to be predicted at. Rows for points, columns for variables.
#' @param fit a list of settings, including the model fit. This can either be a fit created by \code{\link{spotModel.train}} or a spotConfig list as returned by \code{spot}.
#'
#' @return the list of settings, including the model fit. This list can be employed as the "fit" parameter in \code{\link{spotModel.predict}}.
#'
#' @seealso \code{\link{spotModel.train}}, \code{\link{spotModel.func}} 
#'
#' @export
#'
#' @examples
#' ## simple test function
#' sphereFunction <- function(X) rowSums(X^2)
#' ## sample locations
#' Xtrain <- matrix(runif(20),,2)
#' ## evaluate on test function
#' y <- sphereFunction(Xtrain) 
#' ## train model
#' model <- spotModel.train(Xtrain,y,"spotPredictLm")
#' ## evaluate model at a new sample location
#' newy <- spotModel.predict(c(0.5,0.5),model)
###################################################################################
spotModel.predict <- function(X,fit) {
	X <- data.frame(X)
	spotConfig <- eval(call(fit$seq.predictionModel.func, NULL, NULL, X, fit,fit$seq.modelFit))
	spotConfig$seq.largeDesignY # single point prediction
}

###################################################################################
#' Model Prediction Interface
#'  
#' spotMode.func creates a function which represents the chosen surrogate model (see \code{method} parameter).
#' The created function will yield values predicted at given new locations.
#' 
#' @param X Matrix of sample locations. Rows for points, columns for variables.
#' @param y Observations at sample locations (vector)
#' @param method which model function should be used specified as a string, e.g., "spotPredictRandomForest", "spotPredictLm", "spotPredictForrester". Most other "spotPredict*" functions provided by SPOT should work.
#' @param settings additional settings, as used by the corresponding \code{spotPredict*} function. NULL means default model settings are used.
#'
#' @return a function of type y=f(X), where \code{X} is new data, \code{f} is the surrogate model, and \code{y} the predicted values.
#'
#' @seealso Alternatively, \code{\link{spotModel.train}} and \code{\link{spotModel.predict}} can be used to do a similar job.
#'
#' @export
#'
#' @examples
#' ## simple test function
#' sphereFunction <- function(X) rowSums(X^2)
#' ## sample locations
#' Xtrain <- matrix(runif(20),,2)
#' ## evaluate on test function
#' y <- sphereFunction(Xtrain) 
#' ## train model
#' model <- spotModel.func(Xtrain,y,"spotPredictLm")
#' ## evaluate model at a new sample location
#' newy <- model(c(0.5,0.5))
#' ## plot model (based on 10*10 point matrix)
#' #spotSurf3d(model,s=10) #uncomment for testing
###################################################################################
spotModel.func <- function(X,y,method="spotPredictTree",settings=NULL) {
  fit <- spotModel.train(X,y,method=method,settings=settings)
  function(Xnew) spotModel.predict(as.matrix(Xnew), fit)
}
