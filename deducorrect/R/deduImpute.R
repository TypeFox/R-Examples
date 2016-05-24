#' Deductive imputation of numerical or categorical values
#'
#' Based on observed values and edit rules, impute as many variables deductively as possible.
#'
#' @note 
#' When \code{adapt} is not \code{NULL}, values in \code{dat} where \code{adapt==TRUE}
#' are replaced with \code{NA}. The output may therefore contain missings at positions
#' that were previously filled (with wrong values, according to \code{adapt}).
#'
#' @references 
#' T. De Waal, J. Pannekoek and S. Scholtus (2011) Handbook of statistical data editing 
#' Chpt 9.2.1 - 9.2.2
#'
#' @param E An \code{editmatrix} or \code{editarray}
#' @param dat A \code{data.frame}
#' @param adapt (optional) A boolean array of dim(dat), e.g. the result editrules::localizeErrors(E,dat). 
#'      Column names must match those of \code{dat}.
#' @param ... arguments to be passed to \code{\link{solSpace}} (numerical data) or \code{\link{deductiveLevels}} (categorical data)
#'
#' @return A \code{\link{deducorrect-object}}
#'
#' @seealso \code{\link{deductiveZeros}}, \code{\link{solSpace}}, \code{\link{deductiveLevels}}
#'
#' @example ../examples/deduImpute.R
#' @export
deduImpute <- function(E, dat, adapt=NULL, ...){
    if (!is.null(adapt)){
        stopifnot(
            dim(adapt) == dim(dat),
            all(colnames(adapt) %in% names(dat))
        )
    }
    UseMethod('deduImpute')
}

#'
#' If \code{E} is an \code{editset}, imputation based on numerical rules (if any) is performed,
#' and imputations violating extra edits are reverted. Next, this procedure is repeated for
#' pure categorical rules. The results are combined and returned in a \code{deducorrect} object.
#'
#' @method deduImpute editset
#' @rdname deduImpute
#'
#' @export
deduImpute.editset <- function(E, dat, adapt=NULL,...){
    if (!is.null(adapt)){
        N <- colnames(adapt)
        for ( n in N ) dat[adapt[,n],n] <- NA
        adapt <- NULL
    }
    et <- editType(E)
    Em <- NULL
    if ( any(et=='mix'))  Em <- E[et=='mix',]
    NUM <- CAT <- TRUE
    toImpute <- is.na(dat[getVars(E)])
#    if ( !is.null(adapt) ) toImpute <- toImpute | adapt
    toImpute <- rowSums(toImpute)    

    if ( any(et=='num') && !is.null(Em) ){ 
        v1 <- violatedEdits(E,dat)
        dnum <- deduImpute.editmatrix(E$num, dat, adapt, ...)
        v2 <- violatedEdits(E,dat)
        rvt <- apply((!v1 |is.na(v1)) & (v2|is.na(v2)) ,1,any)
        rvt <- rvt[!is.na(rvt)]
        if(any(rvt)) dnum <- revert(dnum,rows=rvt)
    } else if ( any(et=='num') ) {
        dnum <- deduImpute.editmatrix(E$num, dat, adapt,...)
    } else {
        NUM <- FALSE
    }

    icat <- et == 'cat'
    if ( any(icat) && !is.null(Em) ){ 
        v1 <- violatedEdits(E, dat)
        dcat <- deduImpute(reduce(E[icat,]$mixcat), dat, adapt, ...)
        v2 <- violatedEdits(E, dat)
        rvt <- apply((!v1|is.na(v1)) & (v2|is.na(v2)), 1, any)
        rvt <- rvt[!is.na(rvt)]
        if(any(rvt)) dcat <- revert(dcat,rows=rvt)
    } else if ( any(icat) ) {
        dcat <- deduImpute(E$mixcat, dat, adapt, ...)
    } else {
        CAT <- FALSE
    }
        
    if ( !NUM & !CAT) return(newdeducorrect(dat))

    
    if ( NUM & !CAT) return(dnum)

    if ( !NUM & CAT) return(dcat)


    catvar <- as.character(unique(dcat$corrections$variable))
    dat[,catvar] <- as.character(dcat$corrected[,catvar])
    numvar <- as.character(unique(dnum$corrections$variable))
    dat[,numvar] <- dnum$corrected[,numvar]
    
    # this dirty little trick keeps R CMD CHECK from detecting unassigned variables...
    old <- variable <- NULL
    corr <- rbind(
        transform(dnum$corrections,
            variable = as.character(variable)
        ),
        transform(dcat$corrections,
            variable = as.character(variable),
            old=as.character(old),
            new=as.character(new)
        )
    )

    imputations <- dnum$status$imputations + dcat$status$imputations 
    st <- combineStatus(dnum$status$status,dcat$status$status)
    st[st=='partial' & imputations == toImpute & toImpute > 0] <- 'corrected'

    status <- data.frame( 
        status      = st,
        num.status  = dnum$status$status,
        cat.status  = dcat$status$status,
        imputations = imputations
    ) 

    newdeducorrect(
        corrected=dat,
        corrections=corr,
        status=status
    )
}



#' Deductive imputation of categorical data
#'
#' \bold{For categorical data:} The funcion \code{\link{deductiveLevels}} is used to derive
#' deductive imputations for as many fields as possible
#'
#'
#' @method deduImpute editarray
#' @rdname deduImpute
#' @export
deduImpute.editarray <- function(E, dat, adapt=NULL, ...){
   
    if (!is.null(adapt)){
        N <- colnames(adapt)
        for ( n in N ) dat[adapt[,n],n] <- NA
        adapt <- NULL
    }
    vars <- getVars(E)
#    if ( is.null(adapt) ){
        a <- logical(length(vars))
        nCandidates <- rowSums(is.na(dat[,vars,drop=FALSE]))    
#    } else {
#        nCandidates <- rowSums(is.na(dat[,vars,drop=FALSE]) | adapt[,vars,drop=FALSE])
#    }

    nImp <- numeric(nrow(dat))
    X <- t(dat[,vars,drop=FALSE])
    imp <- vector(mode='list',length=nrow(dat))
    for ( i in 1:ncol(X) ){
        x <- X[ ,i]
#        if ( !is.null(adapt) ) a <- adapt[i,vars]
        L <- deductiveLevels(E,x,adapt=a, ...)
        X[names(L),i] <- L
        if ( is.null(L) ){
            imp[[i]] <- character(0)
        } else {
            imp[[i]] <- L
        }
    }
    dat[,vars] <- t(X)

    nImp <- sapply(imp,length)
    # derive deducorrect object
    stat <- status(nrow(dat))
    stat[ nCandidates == 0] <- 'valid'
    stat[ nImp > 0 & nImp < nCandidates] <- 'partial'
    stat[ nImp == nCandidates & nCandidates > 0] <- 'corrected'
    stat[ nImp == 0 & nCandidates > 0 ] <- 'invalid'
    # corrections
    
    xi <- do.call(c,imp)

    corrections <- data.frame(
        row = rep(1:nrow(dat),times=nImp),
        variable = names(xi),
        old = rep(NA,length(xi)), # TODO copy actual old value in case adapt != NULL
        new =  xi)
    newdeducorrect(
        corrected = dat,
        corrections = corrections,
        status = data.frame(status=stat,imputations=nImp)
    )
}






#' Based only equality rules, impute as many values as possible.
#'
#' \bold{For numerical data:} Given (equality) rules and a number of values to impute or adapt, in some cases
#' unique solutions can be derived. This function uses \code{\link{solSpace}} and
#' \code{\link{deductiveZeros}} (iteratively) to determine which values can be imputed
#' deductively. Solutions causing new violations of (in)equality rules are rejected by default by testing
#' if the observed values can lead to a feasible record. This may be switched off by passing
#' \code{checkFeasibility=FALSE}. This may be desirable for performance reasons. If \code{adapt}
#' was computed with an error localization algorithm, such as \code{editrules::localizeErrors}, the 
#' feasibility check is also not nessecary.
#'
#'
#' @method deduImpute editmatrix
#'
#' @param tol tolerance to use in \code{\link{solSpace}} 
#'      and in \code{\link{deductiveZeros}} 
#' @param round should the result be rounded?
#'
#' @rdname deduImpute
#' @export 
deduImpute.editmatrix <- function(E, dat, adapt=NULL, tol=sqrt(.Machine$double.eps), round=TRUE, ...){

    if (!is.null(adapt)){
        N <- colnames(adapt)
        for ( n in N ) dat[adapt[,n],n] <- NA
        adapt <- NULL
    }
    # TODO: change adapt-handling
    vars <- getVars(E)
    X <- t(dat[,vars,drop=FALSE])
    a <- logical(length(vars))
    Xi <- array(NA,dim=c(length(vars),ncol(X)))

    dna <- is.na(dat)
    npre <- rowSums(dna)
    npost <- numeric(nrow(dat))

    for ( i in 1:ncol(X) ){
        x <- X[vars,i]
        nMiss <- sum(is.na(x)) + 1

        while( sum(is.na(x)) < nMiss ){
            nMiss <- sum(is.na(x))
            I <- deductiveZeros(E,x)
            if ( any(I) ) x[I] <- 0
            s <- solSpace(E, x, tol=tol, ...)
            if ( !is.null(s) ){
                u <- rowSums(abs(s$C)) == 0
                x[rownames(s$x0)[u]] <- s$x0[u]
                npost[i] <- npost[i] + sum(u)
            }
        }
        Xi[,i] <- x
    }
    ii <- which( is.na(X[vars,]) & !is.na(Xi))

    corrections <- data.frame(
            row     = (ii-1) %/% nrow(Xi) + 1,   
            variable= vars[(ii-1) %% nrow(Xi) + 1],
            old     = rep(NA,length(ii)),
            new     = Xi[ii]
    )

    X[vars,] <- Xi
    dat[,vars] <- t(X[vars,,drop=FALSE])

    nImp <- npre - npost
    stat <- status(nrow(dat))
    stat[npre  == 0 ]               <- 'valid'
    stat[npost == 0 & npre > 0]     <- 'corrected'
    stat[0 < nImp   & nImp < npre ] <- 'partial'
    stat[npre==npost & npre > 0]    <- 'invalid'

    if (round) dat[,vars] <- round(dat[,vars])
    newdeducorrect(
        corrected   = dat,
        corrections = corrections,    
        status      = data.frame(status=stat, imputations=nImp)
    )
}






