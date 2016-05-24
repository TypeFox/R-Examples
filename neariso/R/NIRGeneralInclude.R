################################################################
###
### This is the main interface function for the FLSA
### it calls the right C++ function depending on whether it is a 2-dimensional
### or a 1-dimensional problem
###
################################################################

neariso = function(y, maxBreaks=100, lambda=NULL)
{
    if(!is.vector(y)) {
        stop("y has to be a vector")
    }
    solObj = .Call("NIR",as.numeric(y), PACKAGE="neariso")
    if(is.null(lambda))
    {
        lambda <- nearisoGetBreakpoints(solObj, maxBreaks)
    }
    else {
        ### find the degrees of freedom for the lambdaVec
    } 
    resExplicit = NIROneDimExplicitSolution(solObj,lambda)

    res <-list(solObj=solObj, lambda=resExplicit$lambda, df=resExplicit$df, beta=resExplicit$beta)
    class(res) <- "nearisoPath" 

    return(res)
}



#################################################################
###
### create a list of the breakpoints
### if there are more than 100, gives 100 evenly spaced over the 
### whole region
###
#################################################################

nearisoGetBreakpoints <- function(nearisoPathObj, maxBreaks=100) {
    if(class(nearisoPathObj) != "nearisoPath" && class(nearisoPathObj) != "nearisoSolObj") {
        stop("nearisoPathObj is not of class nearisoPATH")
    }
   
    solObj <- NULL 
    if(class(nearisoPathObj) == "nearisoPath") {
        solObj <- nearisoPathObj$solObj
    }
    else {
        solObj <- nearisoPathObj
    }

    lambdaList <- sort(unique(solObj$mergeLambda), decreasing=FALSE)
    if(lambdaList[1] == -1) {
        lambdaList[1] <- 0
    }
    else {
        lambdaList <- c(0, lambdaList)
    }

    if(length(lambdaList) >  maxBreaks) {
        selectIndex <- (0:(maxBreaks-1))/(maxBreaks-1) * (length(lambdaList)-1)
        selectIndex <- round(selectIndex)
        lambdaList <- lambdaList[selectIndex + 1]
    }

    return(lambdaList)
}


##################################################################
###
### if an object with the solution tree was returned, this object can
### be used to generate explicit solutions
###
##################################################################

nearisoGetSolution <- function(nearisoPathObj, lambda=nearisoGetBreakpoints(nearisoPathObj))
{
    ### check that lambda is ok
    if(is.null(lambda)) {
        stop("lambda has to be specified")
    }
    lambda = checkLambda(lambda)

    if(class(nearisoPathObj)=="nearisoPath")
    {
        res = NIROneDimExplicitSolution(nearisoPathObj$solObj, lambda)
    }
    else
    {
        stop("nearisoPathObj is not of class nearisoPath")
    }
    nearisoPathObj$beta <- res$beta
    nearisoPathObj$lambda <- res$lambda
    nearisoPathObj$df <- res$df

    return(nearisoPathObj)
}


##################################################################
###
### function that checks if lambda has the right format
###
##################################################################

checkLambda <- function(lambda)
{
    ### check that lambda is a numeric vector, increasing and only non-negative elements
    if(!is.numeric(lambda))
    {
        stop("lambda has to be a numeric vector")
    }
    ### make sure that lambda is non-negative
    if(sum(lambda<0)>0)
    {
        stop("lambda has to be non-negative")
    }
    ### sort lambda, make sure all are unique 
    lambda = sort(unique(lambda))
    return(lambda)
}


####################################################################
###
### get an explicit solution from a FLSA solution object
### this is an internal function and will be hidden in the package
###
####################################################################

NIROneDimExplicitSolution <- function(solObj, lambda)
{
    lambda = checkLambda(lambda)
    lambdaAll <- nearisoGetBreakpoints(solObj, solObj$numVars + 1)
    lambdaAll <- c(lambdaAll, 1e10)
    degFree <- solObj$numVars + 1 - cut(lambda, lambdaAll, right=FALSE, labels=FALSE)

    beta <- t(.Call("NIRexplicitSolution",solObj, lambda, PACKAGE="neariso"))

    return(list(beta=beta, lambda=lambda, df=degFree))
}

