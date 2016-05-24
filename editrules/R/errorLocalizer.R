#' Create a backtracker object for error localization
#' 
#' @section Details:
#' Generate a \code{\link{backtracker}} object for error localization in numerical, categorical, or mixed data.
#' This function generates the workhorse program, called by \code{\link{localizeErrors}} with \code{method=localizer}.
#'
#' The returned \code{\link{backtracker}} can be used to run a branch-and-bound algorithm which finds
#' the least (weighted) number of variables in \code{x} that need to be adapted so that all restrictions 
#' in \code{E} can be satisfied. (Generalized principle of Fellegi and Holt (1976)).
#'
#' The B&B tree is set up so that in in one branche,
#' a variable is assumed correct and its value subsituted in \code{E}, while in the other
#' branche a variable is assumed incorrect and \code{\link[=eliminate]{eliminated}} from \code{E}. 
#' See De Waal (2003), chapter 8 or De Waal, Pannekoek and Scholtus (2011) for 
#' a concise description of the B&B algorithm. 
#'
#' 
#' Every call to \code{<backtracker>$searchNext()} returns one solution \code{list}, consisting of
#' \itemize{
#' \item{w: The solution weight.} 
#' \item{adapt: \code{logical} indicating whether a variable should be adapted (\code{TRUE}) or not}}
#'
#' Every subsequent call leads either to \code{NULL}, in which case either all solutions have been found,
#' or \code{maxduration} was exceeded. The property \code{<backtracker>$maxdurationExceeded} indicates if this is
#' the case. Otherwise, a new solution with a weight \code{w} not higher than the weight of the last found solution
#' is returned.
#' 
#' Alternatively \code{<backtracker>$searchBest()} will return the best solution found within \code{maxduration} seconds.
#' If multiple equivalent solutions are found, a random one is returned.
#'
#' The backtracker is prepared such that missing data in the input record \code{x} is already
#' set to adapt, and missing variables have been eliminated already.
#'
#' The backtracker will crash when \code{E} is an \code{\link{editarray}} and one or more values are
#' not in the datamodel specified by \code{E}. The more user-friendly function \code{\link{localizeErrors}}
#' circumvents this. See also \code{\link{checkDatamodel}}.
#'
#' @section Numerical stability issues:
#' For records with a large numerical range (\emph{eg} 1-1E9), the error locations represent solutions that
#' will allow repairing the record to within roundoff errors. We highly recommend that you round near-zero 
#' values (for example, everything \code{<= sqrt(.Machine$double.eps)}) and scale a record with values larger
#' than or equal to 1E9 with a constant factor.
#'
#' @note This method is potentially very slow for objects of class \code{\link{editset}} that contain 
#'  many conditional restrictions.  Consider using \code{\link{localizeErrors}} with the option 
#'  \code{method="mip"} in such cases.
#'
#'
#'
#' @param E an \code{\link{editmatrix}} or an \code{\link{editarray}}
#' @param x a named numerical \code{vector} or \code{list} (if E is an editmatrix), a named character \code{vector} or \code{list} (if E is an editarray), 
#'   or a named \code{list} if E is an \code{\link[=disjunct]{editlist}} or \code{\link{editset}}.
#'    This is the record for which errors will be localized.
#' @param ... Arguments to be passed to other methods (e.g. reliability weights)
#'
#' @return an object of class \code{\link{backtracker}}. Each execution of \code{$searchNext()} yields a solution
#'      in the form of a \code{list} (see details). Executing \code{$searchBest()} returns the lowest-weight solution.
#'      When multiple solotions with the same weight are found, \code{$searchBest()} picks one at random.
#'
#' @example ../examples/errorLocalizer.R
#' @seealso \code{\link{errorLocalizer_mip}}, \code{\link{localizeErrors}}, \code{\link{checkDatamodel}}, \code{\link{violatedEdits}}, 
#'      
#' @references 
#' I.P. Fellegi and D. Holt (1976). A systematic approach to automatic edit and imputation. 
#' Journal of the American Statistical Association 71, pp 17-25
#'
#' T. De Waal (2003) Processing of unsave and erroneous data.  PhD thesis, Erasmus Research institute 
#' of management, Erasmus university Rotterdam. 
#' http://www.cbs.nl/nl-NL/menu/methoden/onderzoek-methoden/onderzoeksrapporten/proefschriften/2008-proefschrift-de-waal.htm
#' 
#' T. De Waal, Pannekoek, J. and Scholtus, S. (2011) Handbook of Statistical Data Editing. Wiley Handbooks
#' on Survey Methodology.
#'
#' @export
errorLocalizer <- function(E, x, ...){
    UseMethod("errorLocalizer")
}

#'
#' @method errorLocalizer editset
#' @rdname errorLocalizer
#' @export
#'
errorLocalizer.editset <- function(E, x, ...){
    D <- disjunct(E)
    errorLocalizer.editlist(E,x,...)
}


# Localize errors in numerical data
#'
#' @method errorLocalizer editmatrix
#' @param weight a \code{lengt(x)} positive weight vector. The weights are assumed to be in the same order as the variables in \code{x}.
#' @param maxadapt maximum number of variables to adapt
#' @param maxweight maximum weight of solution, if weights are not given, this is equal to the 
#' maximum number of variables to adapt. 
#' @param maxduration maximum time (in seconds), for \code{$searchNext()}, \code{$searchAll()} (not for \code{$searchBest}, use 
#'  \code{$searchBest(maxdration=<duration>)} in stead) 
#' @param tol tolerance passed to \code{link{isObviouslyInfeasible}} (used to check for bound conditions).
#'
#' @rdname errorLocalizer
#' @export
errorLocalizer.editmatrix <- function(
            E, 
            x, 
            weight=rep(1,length(x)), 
            maxadapt=length(x), 
            maxweight=sum(weight),
            maxduration=600,
            tol = sqrt(.Machine$double.eps),
            ...){
    stopifnot(
        is.numeric(weight), 
        all(!is.na(weight)), 
        all(weight>=0), 
        length(weight)==length(x)
    )

    if ( !isNormalized(E) ) E <- normalize(E)
    # missings must be adapted, others still have to be treated.
    adapt <- is.na(x)   
    vars <- getVars(E)
    if (!all(vars %in% names(x)) ) stop('E contains variables not in record')

    o <- order(weight[!adapt], decreasing=TRUE)
    totreat <- names(x)[!adapt][o]

    # only treat variables in occuring in editmatrix.
    totreat <- totreat[totreat %in% vars]
    # if variables do not occur in editmatrix, do not adapt.
    adapt[!(names(adapt) %in% vars)] <- FALSE
    # Eliminate missing variables.
    vs <- names(x)
    elimvars <- vs[vs %in% vars & is.na(x)]
    for (v in elimvars) E <- eliminate.editmatrix(E,v)
    wsol <- min(sum(weight[vars %in% totreat | adapt[vars]]), maxweight)

    cp <- backtracker(
        maxduration=maxduration,
        isSolution = {

            w <- sum(weight[adapt])

            if (.memfail ){
              memfail <- TRUE
              .memfail <- FALSE
              return(FALSE)
            }

            if ( w > min(wsol,maxweight)
              || sum(adapt) > maxadapt
              || isObviouslyInfeasible.editmatrix(.E,tol=tol)
               ) return(FALSE)

            # shortcut: can this ever lead to a solution?
            if ( w == wsol && 
                isObviouslyInfeasible.editmatrix(substValue(.E,totreat,x[totreat]),tol=tol)
            ) return(FALSE)

            if (length(totreat) == 0){
                wsol <<- w
                adapt <- adapt 
                rm(totreat)
                return(TRUE)
            }
        },
        choiceLeft = {
            .var <- totreat[1]
            .E <- substValue.editmatrix(.E, .var , x[[.var]])
            adapt[.var] <- FALSE
            totreat <- totreat[-1]
        },
        choiceRight = {
            .memfail <- FALSE
            .var <- totreat[1]
            .E <- tryCatch(
              eliminate.editmatrix(.E, .var),
              error=function(e){
                .memfail <<- TRUE
                .E 
              })
            adapt[.var] <- TRUE
            totreat <- totreat[-1]
        },
        .E = E,
        x = x,
        maxadapt  = maxadapt,
        maxweight = maxweight,
        totreat   = totreat,
        adapt     = adapt,
        weight    = weight,
        wsol      = wsol,
        tol       = tol,
        .memfail   = FALSE
    )
    
    # add a searchBest function, currently returns random solution in the case of multiple optima
    with(cp,{
        degeneracy <- NA
        searchBest <- function(maxduration=600, VERBOSE=FALSE){
            l <- searchAll(maxduration=maxduration,VERBOSE=VERBOSE)
            if (length(l)>0){ 
                ws <- sapply(l,function(s) s$w)
                iwmin <- which(ws==min(ws))
                degeneracy <<- length(iwmin)
                if (length(iwmin) == 1) return(l[[iwmin]])
                return(l[[sample(iwmin,1)]])
            }
        }
    })
    cp
}


# Localize errors in categorical data
#'
#' @method errorLocalizer editarray
#' @rdname errorLocalizer
#' @export
errorLocalizer.editarray <- function(
    E, 
    x, 
    weight=rep(1,length(x)), 
    maxadapt=length(x),
    maxweight=sum(weight),
    maxduration=600,
    ...){
    
    stopifnot(
        is.numeric(weight), 
        all(!is.na(weight)), 
        all(weight>=0), 
        length(weight)==length(x)
    )
    adapt <- is.na(x)

    vars <- getVars.editarray(E)
    cont <- names(x)[!adapt] %in% vars    
    if (!all(vars %in% names(x)) ) stop('E contains variables not in record')

    o <- order(weight[!adapt], decreasing=TRUE)[cont]
    totreat <- names(x)[!adapt][o]

    for (v in vars[adapt & names(x) %in% vars]) E <- eliminate.editarray(E,v)
    wsol <- min(sum(weight[vars %in% totreat | adapt[vars]]),maxweight)
    ind <- getInd(E)
    bt <- backtracker(
        isSolution = {
            w <- sum(weight[adapt])

            if (.memfail ){
              memfail <- TRUE
              .memfail <- FALSE
              return(FALSE)
            }
            
            if ( w > min(wsol,maxweight) || sum(adapt) > maxadapt )  return(FALSE) 

            # check feasibility 
            .I <- unique(do.call(c, c(ind[.totreat],ind[names(adapt)[adapt]])))
            if ( length(.totreat) > 0 &&  any(rowSums(.E[,.I,drop=FALSE]) == length(.I)) ) return(FALSE)
            
            if ( length(.totreat) == 0 ){
                # can eliminated variables be filled in?
                .I <- do.call(c,ind[adapt])
                if ( length(.I) > 0 && any(rowSums(.E[,.I,drop=FALSE]) == length(.I)) ) return(FALSE)
                # Take care of corner case: check that the record is invalid
                
                if ( length(.I) == 0 && nrow(.E) > 0 && any(apply(.E[,,drop=FALSE],1,all)) )  return(FALSE)
                # prepare output
                wsol <<- w
                adapt <- adapt 
                return(TRUE)
            } 
        },
        choiceLeft = {
            .var <- .totreat[1]
            .E <- substValue.editarray(.E, .var , x[[.var]], reduce=FALSE)
            adapt[.var] <- FALSE
            .totreat <- .totreat[-1]
        },
        choiceRight = {
            .memfail <- FALSE
            .var <- .totreat[1]
            .E <- tryCatch(
              eliminate.editarray(.E, .var),
              error = function(e){
                .memfail <<- TRUE
                .E
              })
            adapt[.var] <- TRUE
            .totreat <- .totreat[-1]
        },
        .E       = E,
        x       = x,
        maxadapt= maxadapt,
        maxweight=maxweight,
        .totreat = totreat,
        adapt   = adapt,
        weight  = weight,
        wsol    = wsol,
        ind     = ind,
        .memfail = FALSE
    )
    # add a searchBest function, currently returns random solution in the case of multiple optima
    with(bt,{
        degeneracy <- NA
        searchBest <- function(maxduration=600, VERBOSE=FALSE){
            l <- searchAll(maxduration=maxduration,VERBOSE=VERBOSE)
            if (length(l)>0){ 
                ws <- sapply(l,function(s) s$w)
                iwmin <- which(ws==min(ws))
                degeneracy <<- length(iwmin)
                if (length(iwmin) == 1) return(l[[iwmin]])
                return(l[[sample(iwmin,1)]])
            }
        }
    })
    bt
}

#'
#' @method errorLocalizer editlist
#' @rdname errorLocalizer
#' @export
errorLocalizer.editlist <- function(
    E, 
    x, 
    weight=rep(1,length(x)), 
    maxadapt=length(x),
    maxweight=sum(weight),
    maxduration=600,
    ...){
    
    stopifnot(
        is.numeric(weight), 
        all(!is.na(weight)), 
        all(weight>=0), 
        length(weight)==length(x) 
    )
    adapt <- is.na(x)

    vars <- getVars.editlist(E)
    cont <- names(x)[!adapt] %in% vars    
    if (!all(vars %in% names(x)) ) stop('E contains variables not in record')

    o <- order(weight[!adapt], decreasing=TRUE)[cont]
    totreat <- names(x)[!adapt][o]

    catvar <- getVars.editlist(E,type='cat')
    icat <- logical(length(adapt))
    names(icat) <- names(adapt)
    icat[catvar] <- TRUE
    for (v in vars[adapt & names(x) %in% vars]) E <- eliminate.editlist(E,v)
    wsol <- min(sum(weight[vars %in% totreat | adapt[vars]]),maxweight)
    # index may be read like this, since the datamodel is constant over all elements of an editlist
    ind <- getInd(E[[1]]$mixcat)
    bt <- backtracker(
        isSolution = {
            w <- sum(weight[adapt])
            if (.memfail ){
              memfail <- TRUE
              .memfail <- FALSE
              return(FALSE)
            }
            
            if ( w > min(wsol,maxweight) || sum(adapt) > maxadapt )  return(FALSE) 

            # check feasibility 
#            .infNum <- sapply(.E,function(e) isObviouslyInfeasible(e$num))
            .inf <- isObviouslyInfeasible(.E)
            if ( all(.inf) ) return(FALSE)
            if ( any(.inf) ) .E <- .E[!.inf]
            # shortcut 1: 
#            if (all(.infNum)) return(FALSE)
            # shortcut 2: bound if numerical part cannot lead to a solution     
            .catvar <- .catvar
            if ( w == wsol ){
                .subvar <- .totreat[!.totreat %in% .catvar]
                if (isObviouslyInfeasible(.E$num, .subvar,x[.subvar])) return(FALSE)
            }


            # categorical variables
            .icat <- .icat
            .infCat <- logical(length(.E))
            .catadapt <- names(adapt)[adapt & .icat]
            if ( length(.catadapt) > 0 ){
                .trcat <- .totreat[.totreat %in% .catvar]
                if ( length(.totreat) > 0 ){
                    .I <- unique( unlist(c(ind[.totreat],ind[.catadapt])) )
                } else  { # nothing to treat...
                    .I <- do.call(c,ind[.catadapt])
                }
                if ( length(.I) > 0 ) {
                    .infCat <- sapply(.E, function(e){
                            any(rowSums(e$mixcat[,.I,drop=FALSE]) == length(.I)) 
                    })
                } #else { # corner case: check that the record is invalid
                  #  .infCat <- sapply(.E, function(e){
                  #      nrow(e$mixcat) > 0 && any(apply(e$mixcat[,,drop=FALSE],1,all)) 
                  #  })
                  #}
            }
            # No feasible region left: bound
            if ( all( .infCat) ) return(FALSE)

            # remove infeasible regions
            .E <- .E[!.infCat]
            
            # at leaf:
            if ( length(.totreat) == 0 ){
                
                # prepare output
                wsol <<- w
                adapt <- adapt 
                return(TRUE)
            } 
        },
        choiceLeft = {
            .var <- .totreat[1]
            .E <- substValue.editlist(.E, .var , x[[.var]], simplify=FALSE, reduce=FALSE)
            adapt[.var] <- FALSE
            .totreat <- .totreat[-1]
        },
        choiceRight = {
            .memfail <- FALSE
            .var <- .totreat[1]
            .E <- tryCatch(
              eliminate.editlist(.E, .var),
              error = function(e){
                .memfail <<- TRUE
                .E
              })
            adapt[.var] <- TRUE
            .totreat <- .totreat[-1]
        },
        .E          = E,
        x           = x,
        maxadapt    = maxadapt,
        maxweight   = maxweight,
        .totreat    = totreat,
        adapt       = adapt,
        weight      = weight,
        wsol        = wsol,
        ind         = ind,
        .icat       = icat,
        .catvar     = catvar,
        .memfail     = FALSE
    )
    # add a searchBest function, returns random solution when multiple weights are encountered.
    with(bt,{
        degeneracy <- NA
        searchBest <- function(maxduration=600, VERBOSE=FALSE){
            l <- searchAll(maxduration=maxduration,VERBOSE=VERBOSE)
            if (length(l)>0){ 
                ws <- sapply(l,function(s) s$w)
                iwmin <- which(ws==min(ws))
                degeneracy <<- length(iwmin)
                if (length(iwmin) == 1) return(l[[iwmin]])
                return(l[[sample(iwmin,1)]])
            }
        }
    })
    bt
}

