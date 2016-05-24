#' Check data against a datamodel
#'
#' Categorical variables in \code{dat} which also occur in \code{E} are checked against the datamodel for 
#' those variables. Numerical variables are checked against edits in \code{E} that contain only a single
#' variable (e.g. \eqn{x > 0}). Values violating such edits as well as empty values are set to adapt.
#'
#' @param E an object of class \code{\link{editset}}, \code{\link{editarray}}, or \code{\link{editmatrix}}
#' @param dat a \code{data.frame}
#' @param weight vector of weigths for every variable of \code{dat} or an array of weight of the same dimensions as \code{dat}.
#' @param ... arguments to be passed to or from other methods
#'  
#' @return An object of class \code{\link{errorLocation}}.
#' @seealso \code{\link{errorLocation}}, \code{\link{localizeErrors}}.
#' @export
checkDatamodel <- function(E, dat,weight=rep(1,ncol(dat)), ...){
    UseMethod('checkDatamodel')
}

# 
#
#' @method checkDatamodel editmatrix
#' @export 
#' @keywords internal
checkDatamodel.editmatrix <- function(E, dat, weight=rep(1,ncol(dat)), call=sys.call(), ...){
    if (nrow(E)==0) return(emptyerrorlocation(dat ,method="checkDatamodel",call=call,... ))
    I <- rowSums( abs(getA(E)) > 1e-8 ) == 1
    localize_singleton(E[I,], dat, weight, method="checkDatamodel", call=call, ...  )
}

# Check categorical data against datamodel
#
#

#' @method checkDatamodel editarray
#' @keywords internal
#' @export 
checkDatamodel.editarray <- function(E, dat, weight=rep(1,ncol(dat)), ...){
    if (any(! (getVars(E) %in% names(dat)) ) ){ 
            vr <- paste(getVars(E)[!getVars(E) %in% names(dat)],collapse=', ')
            stop(paste('Variables',vr,'defined in E not present in dat'))
    }

    m <- nrow(dat)
    if ( is.vector(weight) || (is.array(weight) && nrow(weight)==1) ){
        weight <- t(array(rep(weight,m),dim=c(ncol(dat),m)))
        dimnames(weight) <- dimnames(dat)
    }
    I <- names(dat)[names(dat) %in% getVars(E)]
    adapt <- array(FALSE,dim=dim(dat),dimnames=dimnames(dat))
    ind <- getInd(E)

    w = rep(0,m)
    for ( ii in I ){
        J <- !(dat[,ii] %in% names(ind[[ii]]))
        adapt[,ii] <- J
        w[J] <- w[J] + weight[J,ii]
    }
    status <- emptyStatus(n=m)
    status$weight <- w
    newerrorlocation(
        adapt  = adapt,
        status = status,
        method = 'checkDatamodel',
    )
}


#
# For an \code{\link{editset}}, the categorical variables are tested against
# the categorical datamodel.
#
#

#' @method checkDatamodel editset
#' @export
#' @keywords internal
checkDatamodel.editset <- function(E, dat, weight=rep(1,length(getVars(E))), ...){
    if ( is.null(names(weight)) ) names(weight) <- getVars(E)
    et <- editType(E)
    err <- NULL
    if ( ncol(E$mixcat) > 0 ){
        dvar <- getVars(E,type="dummy")
        iv <- -(1:ncol(E$mixcat))
        if ( length(dvar) > 0 ) iv <- unlist(getInd(E$mixcat)[dvar])
        A <- getArr(E$mixcat)[,-iv,drop=FALSE]
        sep <- getSep(E$mixcat)
        F <- neweditarray(
            A,
            indFromArray(A,sep),
            sep
        )
        err <- checkDatamodel.editarray(F,dat,weight)
    }
    G <- reduce(E$num)
    err %+% checkDatamodel.editmatrix(G,dat,weight)
}





