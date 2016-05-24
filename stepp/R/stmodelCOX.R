#################################################################
#
# stmodelCOX.R
#
#######################
# stepp model: COX    #
#######################
#
setClass("stmodelCOX",
	   representation(MM	   = "ANY"   	# model matrix
				),
	   contains="stmodelKM"
	   )

setMethod("estimate",	
	    signature="stmodelCOX",
	    definition=function(.Object, sp, ...){
		est <- callNextMethod(.Object, sp, ...)

		coltrt    <- .Object@coltrt
		trts      <- .Object@trts
		survTime  <- .Object@survTime
		censor    <- .Object@censor
		covariate <- .Object@MM

		# overwrite KM HR estimate with COX PH HR estimate
		est$model <- "COXe"

		## Create a formula for the Cox model:
		## Cannot handle a no intercept model here
		txassign <- rep(NA, length(coltrt))
 		txassign[which(coltrt == trts[1])] <- 1
 		txassign[which(coltrt == trts[2])] <- 0

		if (is.null(covariate)){
		  fmla <- as.formula(paste("Surv(survTime,censor) ~ txassign"))
		} else 
		if (class(covariate) == "matrix"){
		  if (colnames(covariate)[1]=="(Intercept)"){
		    xnam <- colnames(covariate)[2:dim(covariate)[2]]
		    fmla <- as.formula(paste("Surv(survTime,censor) ~ txassign + ", paste(xnam, collapse= "+")))
		  } else
		  {
		    xnam <- colnames(covariate)[1:dim(covariate)[2]]
		    fmla <- as.formula(paste("Surv(survTime,censor) ~ txassign + ", paste(xnam, collapse= "+")))
		  }
		}
		else
		  stop("Covariates must be specified through a model matrix.")

    		for (i in 1:sp@nsubpop) {
	        Cox	       	<- coxph(fmla, subset=(sp@subpop[,i]==1))
		  if (is.null(covariate)){
	          est$logHR[i]    <- Cox$coefficient
	          est$logHRSE[i]  <- sqrt(Cox$var)
		  } else {
	          est$logHR[i]    <- Cox$coefficient[1]
	          est$logHRSE[i]  <- sqrt(Cox$var[1,1])
		  }
		}
		est$logHRw <- sum(est$logHR/est$logHRSE)
		Cox		       <- coxph(fmla)
		  if (is.null(covariate)){
		    est$overallLogHR   <- Cox$coefficient
		    est$overallLogHRSE <- sqrt(Cox$var)
		  } else {
		    est$overallLogHR   <- Cox$coefficient[1]
		    est$overallLogHRSE <- sqrt(Cox$var[1,1])
		  }

		return(est)
	    }
)

setMethod("test",
	    signature="stmodelCOX",
	    definition=function(.Object, nperm, sp, effect, showstatus=TRUE){

		coltrt    <- .Object@coltrt
		trts      <- .Object@trts
		survTime  <- .Object@survTime
		censor    <- .Object@censor
		covariate <- .Object@MM

		nocovar   <- is.null(covariate)
		result <- callNextMethod(.Object, nperm, sp, effect, showstatus, Cox=TRUE, MM=covariate)
		result$model <- "COXt"
		return(result)
	    }
)

setMethod("print",
	    signature="stmodelCOX",
	    definition=function(x, stobj, estimate=TRUE, cov=TRUE, test=TRUE, ...){
		cat("\n")
      	write(paste("Hazard Ratio Estimates Based on Cox PH Model"), file = "")
		callNextMethod(x, stobj, estimate, cov, test)
	    }
)

# constructor function for stmodelCOX
stepp.COX <- function(coltrt, survTime, censor, trts, timePoint, MM){
	model <- new("stmodelCOX", coltrt=coltrt, survTime=survTime, censor=censor, 
			trts=trts, timePoint=timePoint, MM=MM)
	return(model)
}



