###############################################################
###
### The following function creates a connection list
### for the 2-dimensional Fused Lasso Signal Approximator
###
###############################################################

connListTwoDimensions = function(dimensions)
{
    dimensions=as.integer(dimensions) # need integer
    if(!is.vector(dimensions) || length(dimensions)>2 || (!is.numeric(dimensions) && !is.integer(dimensions)))
    {
        stop("dimensions has to be a numeric vector of length 2")
    }
    
    if(dimensions[1]<2 || dimensions[2]<2)
    {
        stop("Each dimension has to have at least length 2")
    }

    conn = .Call("conn2Dim", dimensions, PACKAGE="PMA")

    ### make the list
    nodeNumbers = matrix(0:(dimensions[1]*dimensions[2]-1), nrow=dimensions[1])
    connList =list(nodes=nodeNumbers, conn=conn)
    class(connList) = "connListObj"
    return(connList)
}


##################################################################
###
### checks an object if it conforms to the specifications of a connection List object
###
##################################################################

is.connListObj = function(obj)
{
    if(class(obj)!="connListObj")
    {
        stop("Object does not have the right class")
    }
    ### check that the nodeNumbers are integers
    if(!is.integer(obj$nodes))
    {
        stop("The nodes part has to be a vector of integers")
    }
    ### check that the conn part has the right length
    if(!(length(obj$conn)==length(obj$nodes)))
    {
        stop("The length of the conn part does not correspond to the number of nodes")
    }
    ### check that every element of the conn part is an integer vector or NULL
    for(i in 1:length(obj$nodes))
    {
        if(!(is.null(obj$conn[[i]]) || is.integer(obj$conn[[i]])))
        {
            stop("All elements of the conn part have to be null or integer vectors")
        }
    }
    ### check that all node numbers that occur in the conn part are eligible
    for(i in 1:length(obj$nodes))
    {
        if(sum(!is.element(obj$conn[[i]], obj$nodes))>0)
        {
            stop(paste("Node",i,"has a connection to a non-existing node"))
        }
    }
    return(TRUE)
}

################################################################
###
### This is the main interface function for the FLSA
### it calls the right C++ function depending on whether it is a 2-dimensional
### or a 1-dimensional problem
###
################################################################

FLSA = function(y, lambda1=0, lambda2=NULL, connListObj = NULL, splitCheckSize=1e9, verbose=FALSE, thr=10e-10, maxGrpNum=4*length(y))
{
    splitCheckSize=as.integer(splitCheckSize)
    if(is.null(connListObj)) # call the appropriate method depending on dimension of y
    {
        if(is.vector(y)) ### call the basic FLSA
        {
            solObj = .Call("FLSA",y, PACKAGE="PMA")
            if(!is.null(lambda2))
            {
                resLambda1Is0 = FLSAOneDimExplicitSolution(solObj,lambda2)
                if(lambda1!=0)
                {
                    res = softThresholding(resLambda1Is0,lambda1)
                    return(res)
                }
                else
                {
                    return(resLambda1Is0)
                }
            }
            else
            {
                return(solObj)
            }
        }
        else if(is.matrix(y)) ### call the 2-dimensional FLSA
        {
            connListObj = connListTwoDimensions(dim(y))
            ### check that everything is ok with lambda2
            if(!is.null(lambda2))
            {
                lambda2 = checkLambda2(lambda2)
                res=.Call("FLSAGeneralMain", connListObj, as.vector(y), lambda2, splitCheckSize, verbose, thr, as.integer(maxGrpNum), PACKAGE="PMA")  
                
                ### format the result in 2 dimensions
                res= array(res, dim=c(length(lambda2), dim(y)))
                ### reset the names of the dimensions
                myDimNames = list(lambda2, 1:(dim(y)[1]), 1:(dim(y)[2]))
                dimnames(res) = myDimNames
                
                ### take lambda1 into account if necessary
                if(lambda1!=0)
                {
                    res = softThresholding(res, lambda1)
                }
            }
            else
            {
                res=.Call("FLSAGeneralMain", connListObj, as.vector(y), lambda2, splitCheckSize, verbose,thr, as.integer(maxGrpNum), PACKAGE="PMA")
            }
            return(res)
        }
    }
    else ### call the general FLSA with the connection list
    {
        if(!is.null(lambda2))
        {
            lambda2=checkLambda2(lambda2)
        }
        if(length(connListObj$nodes)!=length(y))
        {
            stop("y has to have the same number of nodes as connListObj")
        }
        res=.Call("FLSAGeneralMain", connListObj, as.vector(y), lambda2, splitCheckSize, verbose, thr, as.integer(maxGrpNum), PACKAGE="PMA")  
        
        ### take lambda1 into account if necessary
        if(!is.null(lambda2) && (lambda1!=0))
        {
            res = softThresholding(res, lambda1)
        }
        return(res)
    }
}



##################################################################
###
### if an object with the solution tree was returned, this object can
### be used to generate explicit solutions
###
##################################################################

FLSAGetSolution = function(solObj, lambda2, lambda1=0, dim=NULL)
{
    ### check that lambda2 is ok
    lambda2 = checkLambda2(lambda2)

    if(class(solObj)=="FLSA")
    {
        res = FLSAOneDimExplicitSolution(solObj, lambda2)
        if(lambda1!=0)
        {
            res = softThresholding(res, lambda1)
        }
        return(res)
    }
    else if(class(solObj)=="FLSAGeneral")
    {
        ### calculate the explicit solution
        nodes = as.integer(which(solObj$InitialNodeMap>=0)-1)
        res = .Call("FLSAGeneralExplicitSolution",solObj,nodes, lambda2, PACKAGE="PMA")
        
        ###  format in the right way if necessary
        if(!is.null(dim))
        {
            ### check that the dimensions are the same as the number of nodes
            if(prod(dim)!=length(nodes))
            {
                stop("Dimensions are not compatible with solObj")
            }
            res = array(res, dim=c(length(lambda2),dim))
        }
        ### take a look at lambda1
        if(lambda1!=0)
        {
            res = softThresholding(res,lambda1)
        }
        return(res)
    }
    else
    {
        stop("solObj is not of class FLSA or FLSAGeneral")
    }
}


##################################################################
###
### function that checks if lambda2 has the right format
###
##################################################################

checkLambda2 = function(lambda2)
{
    ### check that lambda2 is a numeric vector, increasing and only non-negative elements
    if(!is.numeric(lambda2))
    {
        stop("lambda2 has to be a numeric vector")
    }
    ### make sure that lambda2 is non-negative
    if(sum(lambda2<0)>0)
    {
        stop("lambda2 has to be non-negative")
    }
    ### sort lambda2, make sure all are unique 
    lambda2 = sort(unique(lambda2))
    return(lambda2)
}


####################################################################
###
### get an explicit solution from a FLSA solution object
### this is an internal function and will be hidden in the package
###
####################################################################

FLSAOneDimExplicitSolution = function(solObj, lambda2)
{
    lambda2 = checkLambda2(lambda2)
    return(.Call("FLSAexplicitSolution",solObj, lambda2, PACKAGE="PMA"))
}


###################################################################
###
### Given a matrix with the solutions for lambda2, give back a three-dimensional
### array with soltions for varying lambda1
### This is an internal function 
###
###################################################################

softThresholding = function(solMat, lambda1)
{
    ### check that lambda1 is a numeric vector, increasing and only non-negative elements
    if(!is.numeric(lambda1))
    {
        stop("lambda1 has to be a numeric vector")
    }
    ### make sure that lambda2 is non-negative
    if(sum(lambda1<0)>0)
    {
        stop("lambda1 has to be non-negative")
    }
    ### sort lambda1, make sure all are unique 
    lambda1 = sort(unique(lambda1))

    ### generate the new data array and set dimension and names right
    oldDim = dim(solMat)
    newDim = c(length(lambda1), oldDim)
    oldDimNames = dimnames(solMat)
    if(is.null(oldDimNames))
    {
        oldDimNames = vector("list",2)
    }
    newDimNames = c(list(lambda1), oldDimNames)
    res = array(dim=newDim, dimnames = newDimNames)
    
    ### fill in the soft-thresholded data
    if(length(oldDim)==2)
    {
        for(i in 1:length(lambda1))
        {
            foo = abs(solMat)-lambda1[i]
            foo[foo<0]=0
            res[i,,] = sign(solMat) * foo
        }   
    }
    else if(length(oldDim)==3)
    {
        for(i in 1:length(lambda1))
        {
            foo = abs(solMat)-lambda1[i]
            foo[foo<0]=0
            res[i,,,] = sign(solMat) * foo
        }   
    }
    else
    {
        stop("Wrong dimension of solMat; please inform the maintainer of this package")
    }
    return(res)
}
