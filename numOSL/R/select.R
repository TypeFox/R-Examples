#####
select <-
function(obj, position=NULL, runNumber=NULL, nPoints=NULL, 
         low=NULL, high=NULL,irrTime=NULL, lumType=NULL, 
         dataType=NULL, lightSource=NULL) {
    UseMethod("select")
} #
### 2015.05.01.
select.default <- 
function(obj, position=NULL, runNumber=NULL, nPoints=NULL, 
         low=NULL, high=NULL, irrTime=NULL, lumType=NULL, 
         dataType=NULL, lightSource=NULL) {
    ### Stop if not.
    stopifnot(class(obj)=="binfile", is.list(obj), length(obj)==2L,
              all(names(obj) %in% c("records","tab")),
              length(attributes(obj$records[[1L]]))==9L,
              all(names(attributes(obj$records[[1L]])) %in% 
              c("position", "runNumber", "nPoints","low","high",
              "irrTime","lumType","dataType","lightSource")),
              is.null(position) || (is.numeric(position) && all(position>0L)),
              is.null(runNumber) || (is.numeric(runNumber) && all(runNumber>0L)),
              is.null(nPoints) || (is.numeric(nPoints) && all(nPoints>0L)),
              is.null(low) || (length(low)==1L && is.numeric(low)),
              is.null(high) || (length(high)==1L && is.numeric(high)),
              is.null(irrTime) || (is.numeric(irrTime) && all(irrTime>=0.0)),
              is.null(lumType) || is.character(lumType),
              is.null(dataType) || is.character(dataType),
              is.null(lightSource) || is.character(lightSource))
    ###
    ### Search a position (VECTOR).
    if (!is.null(position))  {
        tab <- obj$tab[(obj$tab[,1L,drop=TRUE] %in% position),,drop=FALSE]
        if (nrow(tab)<=0L)  {
            stop("Error: no record satisfies the given conditions!")
        } # end if.
    } else {
        tab <- obj$tab
    } # end if.
    ###
    ### Search runNumber(s) (VECTOR).
    if (!is.null(runNumber))  {
        tab <- tab[(tab[,2L,drop=TRUE] %in% runNumber),,drop=FALSE]
        if (nrow(tab)<=0L)  {
            stop("Error: no record satisfies the given conditions!")
        } # end if.
    } # end if.
    ###
    ### Search nPoints(s) (VECTOR).
    if (!is.null(nPoints))  {
        tab <- tab[(tab[,3L,drop=TRUE] %in% nPoints),,drop=FALSE]
        if (nrow(tab)<=0L)  {
            stop("Error: no record satisfies the given conditions!")
        } # end if.
    } # end if.
    ###
    ### Search low (DOUBLE).
    if (!is.null(low))  {
        tab <- tab[abs(tab[,4L,drop=TRUE]-low)<=
               sqrt(.Machine$double.eps),,drop=FALSE]
        if (nrow(tab)<=0L)  {
            stop("Error: no record satisfies the given conditions!")
        } # end if.
    } # end if.
    ###
    ### Search high (DOUBLE).
    if (!is.null(high))  {
        tab <- tab[abs(tab[,5L,drop=TRUE]-high)<=
               sqrt(.Machine$double.eps),,drop=FALSE]
        if (nrow(tab)<=0L)  {
            stop("Error: no record satisfies the given conditions!")
        } # end if.
    } # end if.
    ###
    ### Search irrTime(s) (VECTOR).
    if (!is.null(irrTime))  {
        tab <- tab[tab[,6L,drop=TRUE] %in% irrTime,,drop=FALSE]
        if (nrow(tab)<=0L)  {
            stop("Error: no record satisfies the given conditions!")
        } # end if.
    } # end if.
    ###
    ### Search lumType(s) (CHARACTERS).
    if (!is.null(lumType))  {
        tab <- tab[(tab[,7L,drop=TRUE] %in% lumType),,drop=FALSE]
        if (nrow(tab)<=0L)  {
            stop("Error: no record satisfies the given conditions!")
        } # end if.
    } # end if.
    ###
    ### Search dataType(s) (CHARACTERS).
    if (!is.null(dataType))  {
        tab <- tab[(tab[,8L,drop=TRUE] %in% dataType),,drop=FALSE]
        if (nrow(tab)<=0L)  {
            stop("Error: no record satisfies the given conditions!")
        } # end if.
    } # end if.
    ###
    ### Search lightSource(s) (CHARACTERS).
    if (!is.null(lightSource))  {
        tab <- tab[(tab[,9L,drop=TRUE] %in% lightSource),,drop=FALSE]
        if (nrow(tab)<=0L)  {
            stop("Error: no record satisfies the given conditions!")
        } # end if.
    } # end if.
    ###
    ###
    recordIndex <- as.numeric(rownames(tab))
    if (all(diff(tab[,3L,drop=TRUE])==0L))  {
        records <- c()
        for (i in recordIndex)  {
            records <- cbind(records, as.numeric(obj$records[[i]]))
        } # end for.
    } else {
        records <- obj$records[recordIndex]
    } # end if.
    ###
    ###
    output <- list("records"=records, "tab"=tab)
    return(output)
} # end function select.
#####
