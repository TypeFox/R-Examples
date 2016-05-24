##################################################
###
### provides the interface functions in R to access
### the C++ functions. 2 are available, one for the 
### whole covariance matrix, one for the covariance
### for just one variable
### Variables are numbered 1 to n
### Also, the log of the partition function can be computed
###
##################################################

BMNJT = function(thetaMat, adjMat=NULL, var=NULL, onlyActive=FALSE, timeout=60)
{
    timeout = as.integer(timeout)
    ### check the input format
    if(!is.matrix(thetaMat))
    {
        stop("thetaMat has to be a matrix")
    }
    if(!(dim(thetaMat)[1]==dim(thetaMat)[2]))
    {
        stop("thetaMat has to be a square matrix")
    }
    if(sum(abs(thetaMat - t(thetaMat)))/length(thetaMat)>10^(-6))
    {
        stop("Theta has to be symmetric")
    }

    if(is.null(adjMat))
    {
        adjMat = (thetaMat!=0)
    }
    
    if(is.null(var)) ### produce complete covariance matrix
    {
        if(onlyActive)
        {
            res = .Call("runJTAlgSecMomMatActive", adjMat, thetaMat, timeout, PACKAGE="BMN")
            return(list(Expectation=diag(res), SecondMomentMatrix=res))
        }
        else
        {
            return(.Call("runJTAlgSecMomMat", adjMat, thetaMat, timeout, PACKAGE="BMN"))
        }
    }
    else
    {
        if(onlyActive)
        {
            stop("Cannot pick a specific variable and onlyActive")
        }
        # also has to be of size 1
        if(length(var)>1)
        {
            stop("var is only allowed to be a single variable")
        }
        size = dim(thetaMat)[1]
        if(var<1 || var>size)
        {
            stop("variable numbers have to be between 1 and the dimension of thetaMat")
        }
        ### internally, variables are numbered 0 to n-1
        var = as.integer(var-1)
        return(.Call("runJTAlgSecMomVec", adjMat, thetaMat, var, timeout, PACKAGE="BMN"))
    }
}

BMNJTlogPartFunc=function(thetaMat, timeout=1e9)
{
    ### check the input format
    if(!is.matrix(thetaMat))
    {
        stop("thetaMat has to be a matrix")
    }
    if(!(dim(thetaMat)[1]==dim(thetaMat)[2]))
    {
        stop("thetaMat has to be a square matrix")
    }
    if(sum(abs(thetaMat - t(thetaMat)))/length(thetaMat)>10^(-6))
    {
        stop("Theta has to be symmetric")
    }

    adjMat = (thetaMat!=0)
  
    return(log(.Call("runJTAlgNormalizationConstant", adjMat, thetaMat, as.integer(timeout), PACKAGE="BMN")))
}

