##Why this file?
##This set of functions is meant to better support DS Calculus
##DS is how we store the "knowledge" that we gain from the data
##The base version of DS was DS_Calculus3a.R, specifically the function (class) bpa
##This base class is very well written, and exposes the basic kinds of operations that we need to support
##However, a single bpa object is capable of representing only the beliefs about a single point
##The classifier is so written that it can deal with a batch of points at the same time
##So we would like the underlying DS objects to be capable of the same

##Features that we want
##Support a batch of points at the same time - stored as a matrix
##Each point is a column, that is each column represents a single bpa object
##Further, with each such matrix, we also have an 'info' field
##This field is used to store the (random seed of the) matrix used to obtain the current bpa's
##If a single point is queried for, as it usually will be, then that column is returned as a bpa object
##Finally, matrices can be combined, as can any DS object

bpamat <- function(info = NULL, mat = NULL)
    {
        info <- info
        mat <- mat
        setlist <- NULL
        points <- NULL

        assign.mat <- function(mat)
            {
                mat <<- mat
                setlist <<- row.names(mat)
                points <<-  colnames(mat)
            }
        
        if(!is.null(mat))
            assign.mat(mat)
        
        get.classify <- function()
            {
                predicted <- apply(mat, MARGIN = 2, FUN = function(x) names(which.max(x)))
                names(predicted) <- points
                return(predicted)
            }

        get.point <- function(p)
            {
                y <- mat[,p]
                b <- bpa(n = dim(mat)[1], setlist = names(y), mlist = unname(y))
                b$set.name(p)
                b$set.info(info)
                return(b)
            }

        ret <- list(set.info = function(info) info <<- info,
                    get.info = function() info,
                    assign.mat = assign.mat,
                    get.classify = get.classify,
                    get.point = get.point,
                    get.mat = function() mat,
                    get.setlist = function() setlist,
                    get.pointlist = function() points)
        class(ret) <- 'bpamat'
        return(ret)
    }

