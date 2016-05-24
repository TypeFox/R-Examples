
#' Check data against constraints
#'
#' Determine which record violates which edits. Returns \code{NA} when edits
#' cannot be checked because of missing values in the data. 
#'
#' @example ../examples/violatedEdits.R
#' @export
#' @seealso \code{\link{checkDatamodel}}
#' @param E \code{\link{character}} vector with constraintsm, \code{\link{editset}}, \code{\link{editmatrix}} or \code{\link{editarray}}.
#' @param dat \code{data.frame} with data that should be checked, if a named vector is supplied it will converted internally to a \code{data.frame}
#' @param ... further arguments that can be used by methods implementing this generic function
#' @return An object of class \code{violatedEdits}, which is a logical \code{nrow(dat)Xnedits(E)} matrix with an extra \code{class} attribute
#'  for overloading purposes. 
violatedEdits <- function(E, dat, ...){
    if (nedits(E)==0){
        v <- matrix( logical(0)
                   , nrow=ifelse(is.vector(dat),1,nrow(dat))
                   , ncol=0
                   )
        dimnames(v) <- list(record=rownames(dat),edit=NULL)
        return(newviolatedEdits(v))
    }
    UseMethod("violatedEdits")
}


#' @rdname violatedEdits
#' @method violatedEdits character
#' @param name name of edits
#' @export
violatedEdits.character <- function(E, dat, name=NULL, ...){
    ed <- parseEdits(E)
    if (is.vector(dat) && !is.null(names(dat))){
       dat <- data.frame(t(dat))
    }
    M <- tryCatch(sapply(ed, eval, envir=dat), error=function(e){
        stop(paste("Not all edits can be evaluated, parser returned", e$message, sep="\n"))})
    if ( is.vector(M) )  M <- array(M, dim=c(1,length(M)))
    dimnames(M) <- list(record=rownames(dat),edit=names(E))
    return(newviolatedEdits(!M))
}

#' Method for editmatrix
#'
#' \itemize{
#' \item{For rules of the form Ax == b  |Ax - b| <= tol is returned.}
#' \item{For rules of the form Ax < b, Ax - b < tol is returned.}
#' \item{For rules of the form Ax <= b Ax- b <= tol is returned.}
#'}
#' For numerical records, the default tolerance is 0. When working with doubles, 
#' the square root of machina accuracy is a resonable alternative (\code{sqrt(.Machine\$double.eps)}).
#' The editmatrix is \cite{\link[=normalize]{normalized}} before checks are performed.
#'
#' @rdname violatedEdits
#' @method violatedEdits editmatrix
#' @param tol tolerance to check rules against.
#' @export
violatedEdits.editmatrix <- function(E, dat, tol=0, ...){
    if (tol < 0 ) stop("Argument tol must be nonnegative")
    if (tol==0) return(violatedEdits.character(as.character(E),dat))
    if ( !isNormalized(E) ) E <- normalize(E)

    eq <- getOps(E) == '=='
    iq <- !eq
    v <- matrix(FALSE,
        nrow=ifelse(is.vector(dat),1,nrow(dat)),
        ncol=nrow(E))

    if ( any(iq) ){
        E[iq,ncol(E)] <- E[iq,ncol(E)] + tol 
        v[,iq] <- violatedEdits.character(as.character(E[iq,]),dat)[,,drop=FALSE]
    }
    if ( any(eq) ){
        nc <- ncol(E)    

        Emin <- E[eq,]
        Emin[,nc] <- Emin[,nc] - tol
        emin <- gsub('==','>=',as.character(Emin))

        Emax <- E[eq,]
        Emax[,nc] <- Emax[,nc] + tol    
        emax <- gsub('==','<=', as.character(Emax))
        v[,eq] <- violatedEdits.character(emin,dat) | violatedEdits.character(emax,dat)
    }

    dimnames(v) <- list(record=rownames(dat),edit=rownames(E))
 
    newviolatedEdits(v)
}

## TODO: Remove this function (mvdl):

# @rdname violatedEdits
# @method violatedEdits data.frame
# @export
violatedEdits.data.frame <- function(E, dat, ...){
    if ( !all(c("name","edit","description") %in% names(E)) ){
        stop("Invalid input data.frame see ?editmatrix for valid input format")
    }
    return(violatedEdits.character(as.character(E$edit), dat, E$name))
}

#'
#'
#' @method violatedEdits editarray
#' @param datamodel Also check against datamodel?
#' @rdname violatedEdits
#' @export
violatedEdits.editarray <- function(E, dat, datamodel=TRUE,...){
    edits <- as.character(E, useIf=FALSE, datamodel=datamodel)
    violatedEdits.character(edits,dat,...)
}

#'
#' 
#' @method violatedEdits editset
#' @rdname violatedEdits
#' @export
violatedEdits.editset <- function(E, dat, datamodel=TRUE, ...){
    v1 <- violatedEdits(E$num,dat,...)
    E$num <- editmatrix(expression())
    u <- as.character(E,datamodel=datamodel,useIf=FALSE)
    v2 <- violatedEdits(u,dat)
    newviolatedEdits(cbind(v1,v2))
}


# TODO: remove the function below. (mvdl)

# Lists which rows of \code{data.frame dat} violate which constraints
#
# This function can be used as an input for automatic corrections methods.
# @example ../examples/listViolatedEdits.R
# 
# @param E a number of edit restrictions, represented as \code{character} vector, \code{\link{editmatrix}} or \code{data.frame}.
# @param dat \code{data.frame} with data that should be checked
# 
# @return a list with per row a \code{integer} vector of the constraints that are violated 
listViolatedEdits <- function(E, dat){    
    errors <- violatedEdits(E, dat)
    errorlist <- apply(errors, 1, which)
    return(apply(errors, 1, which))
}

newviolatedEdits <- function(x){
  dimnames(x) <- list(record=rownames(x),edit=colnames(x))
  structure(x, class="violatedEdits")
}

#' Plot summary statistics on violatedEdits
#' @method plot violatedEdits
#' @param x \code{violatedEdits} object.
#' @param topn Top \code{n} edits to be plotted.
#' 
#' @rdname violatedEdits
#' @export
plot.violatedEdits <- function(x, topn=min(10,ncol(x)), ...){
  N <- nrow(x)
  Nna <- sum(apply(is.na(x),1, all))

  editfreq <- sort(colSums(x, na.rm=TRUE), decreasing=TRUE)
  editfreq <- editfreq[1:topn]
  editfreq <- editfreq/N

  xlim <- c(0, max(editfreq))

  oldpar <- par(mfrow=c(2,1))
  barplot( sort(editfreq, decreasing=TRUE), 
         , main=paste("Edit violation frequency of top",topn,"edits")
         , xlab = "Frequency"
         , ylab= "Edit"
         , horiz = TRUE
         , las = 1
         , xlim = xlim
         , ...
         )
  
  x[is.na(x)] <- TRUE
  cnt <- table(rowSums(x))
  ner <- as.integer(names(cnt))
  cnt <- as.integer(cnt)
  
  noerr <- ner==0
  nnoer <- sum(cnt[noerr],0)
  ner <- ner[!noerr]
  cnt <- cnt[!noerr]
  lgcrit <- 50
  lg <- ''
  if ( max(ner) > lgcrit ) lg <- paste(lg,'x',sep='')
  if ( max(cnt) > lgcrit ) lg <- paste(lg,'y',sep='')
  plot( ner,cnt,
         , main=  "Edit violations per record"
         , xlab = "Number of violations"
         , ylab = "Count"
         , log=lg
         , ...
         )
  mtext(paste(nnoer,'records with no violations'),side=3,line=0)
  par(oldpar)
}

#'
#' @method summary violatedEdits
#' @param object \code{violatedEdits} object
#' @param minfreq minimum freq for edit to be printed
#'
#' @note When summarizing an object of class \code{violatedEdits}, every 
#' empty value is counted as one edit violation when counting violations per record.
#'
#' @rdname violatedEdits
#' @export
summary.violatedEdits <- function(object, E=NULL, minfreq=1, ...){
  N <- nrow(object)
  if (sum(object, na.rm=TRUE)==0){
    message(sprintf("No violations detected, %d checks evaluated to NA",sum(is.na(object))))
    return()
  }
  Nna <- sum(apply(is.na(object),1, all))
  
  editfreq <- sort(colSums(object, na.rm=TRUE), decreasing=TRUE)
  editperc <- round(100*editfreq/N, 1)
  
  editperc <- editperc[editfreq >= minfreq]
  editfreq <- editfreq[editfreq >= minfreq]
  
  editdf <- data.frame(editname=names(editfreq), freq=editfreq, rel=paste(editperc,"%", sep=""))
  if (!is.null(E)){
     editdf$edit <- as.character(E)[editdf$editname]
  }
  
  cat("Edit violations, ",N," observations, ",Nna," completely missing (",round(100*Nna/N,1),"%):\n\n", sep="")
  print(editdf, row.names=FALSE)
  
  object[is.na(object)] <- TRUE
  errfreq <- unclass(table(rowSums(object)))
  errperc <- round(100*errfreq/(N), 1)
  errdf <- data.frame(errors=names(errfreq), freq=errfreq, rel=paste(errperc,"%", sep=""))
  cat("\nEdit violations per record:\n\n")
  print(errdf, row.names=FALSE)
}

#' Print violatedEdits
#' @method print violatedEdits
#' @param x violatedEdits
#' @export
#' @keywords internal
print.violatedEdits <- function(x,...){
  print(unclass(x))
}

#' as.data.frame violatedEdits
#' @method as.data.frame violatedEdits
#' @rdname violatedEdits 
#' @export
as.data.frame.violatedEdits <- function(x, ...){
  as.data.frame(unclass(x),...)
}
