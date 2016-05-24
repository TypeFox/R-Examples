### Constructor (Accepts either controlPoints or parameters)
"SimilarityTransformation" <- function(controlPoints=data.frame(), parameters=numeric()) {
    if (missing(parameters))
        return(new("SimilarityTransformation", controlPoints=controlPoints))
    if (missing(controlPoints))
        return(new("SimilarityTransformation", parameters=parameters))
    new("SimilarityTransformation", controlPoints=controlPoints, parameters=parameters)
}

### Calculate parameters from control points
###     Modifies the original object
### Arguments:
### - object is an "SimilarityTransformation" object
###
if (!isGeneric("calculateParameters"))
    setGeneric("calculateParameters",function(object){standardGeneric ("calculateParameters")})
setMethod(f="calculateParameters",signature(object="SimilarityTransformation"),
    definition=function(object){
		if (ncol(object@controlPoints) == 0)
			stop("Control points were not provided. You could access the parameters directly by calling 'getParameters'.")

        newObject <- deparse(substitute(object))

        x1=c(object@controlPoints[,1],object@controlPoints[,2])
        x2=c(object@controlPoints[,2],-object@controlPoints[,1])
        ones = rep(1,nrow(object@controlPoints))
        zeros = rep(0,nrow(object@controlPoints))
        x3=c(ones,zeros)
        x4=c(zeros,ones)
        linMod=lm(formula=c(object@controlPoints[,3],object@controlPoints[,4]) ~ x1 + x2 + x3 + x4 - 1)

        object@parameters <- as.vector(coef(linMod))
        object@residuals <- matrix(linMod$residuals,nrow(object@controlPoints))
        object@rmse <- sqrt(sum(object@residuals**2)/nrow(object@residuals))

        assign(newObject,object,envir=parent.frame())
		return(invisible())
	}
)

