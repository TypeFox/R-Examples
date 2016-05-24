#' @rdname ForecastData
setGeneric("setPredCalibration<-",function(object,value){standardGeneric("setPredCalibration<-")})

#' @rdname ForecastData
#' @export
setReplaceMethod(
	f="setPredCalibration",
	signature="ForecastData",
	definition=function(object,value){
          if(class(value)=="data.frame"){value <- as.matrix(value)}
          if(class(value)=="matrix"){value <- array(value, dim=c(nrow(value), ncol(value), 1))}
          object@predCalibration = value
          validObject(object)
          return(object)
	}
)


#' @rdname ForecastData
#' @export
setGeneric("setPredTest<-",function(object,value){standardGeneric("setPredTest<-")})

#' @rdname ForecastData
#' @export
setReplaceMethod(
	f="setPredTest",
	signature="ForecastData",
	definition=function(object,value){
          if(class(value)=="data.frame"){value <- as.matrix(value)}
          if(class(value)=="matrix"){value <- array(value, dim=c(nrow(value), ncol(value), 1))}
          object@predTest<- value
		validObject(object)
		return(object)
	}
)

#' @rdname ForecastData
#' @export
setGeneric("setOutcomeCalibration<-",function(object,value){standardGeneric("setOutcomeCalibration<-")})

#' @rdname ForecastData
#' @export
setReplaceMethod(
	f="setOutcomeCalibration",
	signature="ForecastData",
	definition=function(object,value){
		object@outcomeCalibration <- value
		validObject(object)
		return(object)
	}
)

#' @rdname ForecastData
#' @export
setGeneric("setOutcomeTest<-",function(object,value){standardGeneric("setOutcomeTest<-")})

#' @rdname ForecastData
#' @export
setReplaceMethod(
	f="setOutcomeTest",
	signature="ForecastData",
	definition=function(object,value){
		object@outcomeTest<-value
		validObject(object)
		return(object)
	}
)

#' @rdname ForecastData
#' @export
setGeneric("setModelNames<-",function(object,value){standardGeneric("setModelNames<-")})


#' @rdname ForecastData
#' @export
setReplaceMethod(
	f="setModelNames",
	signature="ForecastData",
	definition=function(object,value){
		object@modelNames <-value
		colnames(object@predCalibration)<-value
		colnames(object@predTest)<-value
		validObject(object)
		return(object)
	}
)




