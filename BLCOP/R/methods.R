################ 
#Methods
################

show.BLViews <- function(object)
{
    .innerViewString <- function(p) {
        # p should be a row in the "P" matrix
        # only nonzero elements really consitute "views"
        qv <- tail(p,1)
        x <- head(p,-1)[head(p,-1) != 0]
        tmp <- paste(as.character(x), names(x), sep = "*", collapse = "+")
        tmp <- paste(tmp, qv, sep = "=")
        
    }
    viewStrings <- apply(viewMatrix(object),1,.innerViewString)
    for(i in seq(along  = viewStrings))
        cat(i, ":", viewStrings[i], " + eps. Confidence:", object@confidences[i], " \n")
}

show.BLResult <- function(object)
{
    cat("Prior means:\n")
    show(object@priorMean)
    cat("Posterior means:\n")
    show(object@posteriorMean)
    cat("Posterior covariance:\n")
    show(object@posteriorCovar)
}

show.COPResult <- function(object)
{

	cat(paste("Asset set: ", paste(assetSet(priorViews(object)), collapse = ","), "\n"))
	cat("Views used to generate this posterior: \n")
	show(priorViews(object))
	cat("Number of simulations:", numSimulations(object), "\n" )
	
}

show.COPViews <- function(object)
{
    for(i in 1:nrow(object@pick))
    {                
        x <- object@pick[i, object@pick[i,] != 0]        
        
        tmp <- paste("(", paste(names(object@viewDist[[i]]@parameters), object@viewDist[[i]]@parameters, collapse = ",", sep = "="), ")", sep="")
        distString <- paste(object@viewDist[[i]]@RName, tmp, sep = ":")
        tmp <- paste(as.character(x), names(x), sep = "*", collapse = "+")
        print(paste(tmp, distString, sep = "~"))        
    }    
}

setMethod("show", signature(object = "BLViews"), show.BLViews)
setMethod("show", signature(object = "BLResult"), show.BLResult)
setMethod("show", signature(object = "COPViews"), show.COPViews)
setMethod("show", signature(object = "COPResult"), show.COPResult)