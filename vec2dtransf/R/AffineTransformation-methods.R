### Constructor (Accepts either controlPoints or parameters)
"AffineTransformation" <- function(controlPoints=data.frame(), parameters=numeric()) {
    if (missing(parameters))
        return(new("AffineTransformation", controlPoints=controlPoints))
    if (missing(controlPoints))
        return(new("AffineTransformation", parameters=parameters))
    new("AffineTransformation", controlPoints=controlPoints, parameters=parameters)
}

### Calculate parameters from control points
###     Modifies the original object
### Arguments:
### - object is an "AffineTransformation" object
###
if (!isGeneric("calculateParameters"))
    setGeneric("calculateParameters",function(object){standardGeneric ("calculateParameters")})
setMethod(f="calculateParameters",signature("AffineTransformation"),
    definition=function(object){
		if (ncol(object@controlPoints) == 0)
			stop("Control points were not provided. You could access the parameters directly by calling 'getParameters'.")

        newObject <- deparse(substitute(object))
        names(object@controlPoints) <- c('X_Source','Y_Source','X_Target','Y_Target')   
        linMod <- lm(formula = cbind(X_Target, Y_Target) ~ X_Source + Y_Source, 
                       data=object@controlPoints)

        object@parameters <- as.vector(rbind(coef(linMod)[2:3,],coef(linMod)[1,]))
        object@residuals <- as.matrix(linMod$residuals)
        object@rmse <- sqrt(sum(object@residuals**2)/nrow(object@residuals))

        assign(newObject,object,envir=parent.frame())
		return(invisible())
	}
)

