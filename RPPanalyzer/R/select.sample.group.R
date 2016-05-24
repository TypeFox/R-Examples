#'
#' The function selects subgroups of samples based on different parameters which are listed in the sampledescription
#' @param params: a list of sampledescription column names and associated values
#' @combine: logical value indicating if the union (combine is TRUE) or the intersect (combine is FALSE) should be considered. The default value is FALSE.
#'
`select.sample.group` <-
function(x,params=list(tissue=c("T", "N")), combine=F) {

    # check if every params are column names of the sampledescription
    if(!all(names(params) %in% colnames(x[[4]]))) {

        stop("Following parameters are not contained in the sampledescription: ", paste(names(params)[!(names(params) %in% colnames(x[[4]]))],collapse=" "), "\n")
    
    }

    # create an initial index vector depending if the 
    # the boolean indices should be combined (logical OR) or intersect (logical AND)
    if(combine) {
        dat.lines <- rep(F, nrow(x[[4]]))
    }
    else {
        dat.lines <- rep(T, nrow(x[[4]]))
    }

    # iterate over all given column names
    for (p in names(params)) {

        temp.lines <- x[[4]][,p] %in% params[[p]]  

        if (combine) {
            dat.lines <- dat.lines | temp.lines
        }
        else {
            dat.lines <- dat.lines & temp.lines
        }
    }

    # use th index vector to filter the matrix with the expression values and the variances
    x[[1]] <- x[[1]][dat.lines,]
    x[[2]] <- x[[2]][dat.lines,]
    x[[4]] <- x[[4]][dat.lines,]

    # return the filtered matrix
    return(x)

}

