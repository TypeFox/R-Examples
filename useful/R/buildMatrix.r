#' @title build.x
#' 
#' @description Build the x matrix for a glmnet model
#' 
#' @details Given a formula and a data.frame build the predictor matrix
#' @author Jared P. Lander
#' @aliases build.x
#' @export build.x
#' @importFrom plyr catcolwise
#' @importFrom stats model.matrix terms
#' @param formula A formula
#' @param data A data.frame
#' @param contrasts Logical indicating whether a factor's base level is removed.  Can be either one single value applied to every factor or a value for each factor.  Values will be recycled if necessary.
#' @param sparse Logical indicating if result should be sparse.  Currently not used.
#' @return A matrix of the predictor variables specified in the formula
#' @examples
#' require(ggplot2)
#' head(mpg)
#' head(build.x(hwy ~ class + cyl + year, data=mpg))
#' 
#' testFrame <- data.frame(First=sample(1:10, 20, replace=TRUE), 
#' Second=sample(1:20, 20, replace=TRUE), 
#' Third=sample(1:10, 20, replace=TRUE), 
#' Fourth=factor(rep(c("Alice","Bob","Charlie","David"), 5)), 
#' Fifth=ordered(rep(c("Edward","Frank","Georgia","Hank","Isaac"), 4)), 
#' Sixth=factor(rep(c("a", "b"), 10)), stringsAsFactors=F)
#' head(build.x(First ~ Second + Fourth + Sixth, testFrame, 
#' contrasts=c("Fourth"=TRUE, "Fifth"=FALSE, "Sixth"=TRUE)))
#' head(build.x(First ~ Second + Fourth + Fifth + Sixth, testFrame, 
#' contrasts=c(Fourth=TRUE, Fifth=FALSE, Sixth=TRUE)))
#' head(build.x(First ~ Second + Fourth + Fifth + Sixth, testFrame, contrasts=TRUE))
#' head(build.x(First ~ Second + Fourth + Fifth + Sixth, testFrame, 
#' contrasts=FALSE))
#' head(build.x(First ~ Second + Fourth + Fifth + Sixth - 1, testFrame, 
#' contrasts=TRUE))
#' head(build.x(First ~ Second + Fourth + Fifth + Fourth*Sixth, testFrame, contrasts=TRUE))
#' head(build.x(First ~ Second + Fourth + Fifth + Third*Sixth, testFrame, contrasts=TRUE))
#' #' head(build.x(First ~ Second + Fourth + Fifth + Fourth*Sixth, testFrame, contrasts=FALSE))
#' head(build.x(First ~ Second + Fourth + Fifth + Third*Sixth, testFrame, contrasts=FALSE))
#' 
#' ## if contrasts is a list then you can specify just certain factors
build.x <- function(formula, data, contrasts=TRUE, sparse=FALSE)
{
    # ensure data is a data.frame
    data <- ForceDataFrame(data)
    
    if(length(contrasts) == 1 && contrasts)
    {
        return(model.matrix(formula, data=data))
    }
        
    # make index of factor or character columns
    catIndex <- which(sapply(data, function(x) is.factor(x) | is.character(x)))
    # only keep those that also appear in the factors attr of the terms of formula
    theTerms <- rownames(attr(terms(formula, data=data), "factors"))
    # new cat index only keeping those variables that are necessary
    catIndex <- catIndex[which(names(data)[catIndex] %in% theTerms)]
    # also cut down contrasts argument
    # save for another time
    
    if(length(catIndex) == 0)
    {
        return(model.matrix(formula, data=data))
    }
    
    # if any of these identified columns is still a character, they need to be changed into a factor
    # find out which columns are characters
    #print(sapply(data[, catIndex], is.character))
    charIndex <- catIndex[sapply(data[, catIndex], is.character)]
    
    # convert to factor
    data[, charIndex] <- catcolwise(as.factor)(data[, charIndex, drop=FALSE])
    ## now all factors or characters are at least factors (and nothing extraneous was done) and only the appropriate columns will be put into the contrasts argument
    
    # if multiple contrasts are given they must be named
    contrNames <- names(contrasts)
    if(length(contrasts) > 1 && is.null(contrNames))
    {
        stop("If specifying more than one contrasts then it must be a named list of vector.")
    }else if(!is.null(contrNames))
    {
        # get names of contrasts and use as the catIndex, factor/ordered columns not specified will be left to the default
        catIndex <- contrNames
    }else if(length(contrasts) == 1)
    {
        # make as many contrasts as necessary
        contrasts <- rep(contrasts, times=length(catIndex))
    }
    
    # only non sparse is allowed for now
    sparse <- FALSE
    # build contrast argument list
    #contrArgs <- lapply(data[, catIndex, drop=FALSE], contrasts, contrasts=contrasts, sparse=sparse)
    contrArgs <- mapply(contrasts, data[, catIndex, drop=F], contrasts, MoreArgs=list(sparse=sparse))
    
    # build model.matrix
    model.matrix(formula, data=data, contrasts.arg=contrArgs)#[, -1])
    #model.matrix(formula, data=data)[, -1]
}
#mapply(function(input, contrasts, sparse=FALSE){ contrasts(x=input, contrasts=contrasts, sparse=sparse) }, testFrame[, 4:5, drop=F], c(T), MoreArgs=list(sparse=F))
#head(model.matrix(~ ., data=testFrame, ))
#' ForceDataFrame
#' 
#' Force matrix and arrays to data.frame
#' 
#' This is a helper function for build.x and build.y to convert arrays and matrices--which are not accepted in model.frame--into data.frames
#' 
#' @author Jared P. Lander
#' @aliases ForceDataFrame
#' @return a data.frame of the data
#' @param data matrix, data.frame, array, list, etc. . .
#' 
ForceDataFrame <- function(data)
{
    if(class(data) %in% c("matrix", "array"))
    {
        return(as.data.frame(data))
    }
    return(data)
}

#' build.y
#' 
#' Build the y survival object for a glmnet model
#' 
#' Given a formula and a data.frame build the y survival object
#' @author Jared P. Lander
#' @aliases build.y
#' @export build.y
#' @importFrom stats model.frame
#' @param formula A formula
#' @param data A data.frame
#' @return A surival object for the portion of the formula in Surv
#' @examples
#' require(ggplot2)
#' head(mpg)
#' head(build.y(hwy ~ class + cyl + year, data=mpg))
#' 
build.y <- function(formula, data)
{
    # build a model frame
    theFrame <- model.frame(formula, data=ForceDataFrame(data))
    # extract the response
    theFrame[[1]]
}
