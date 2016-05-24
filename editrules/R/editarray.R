#' Parse textual, categorical edit rules to an editarray
#'
#' An \code{editarray} is a boolean array (with some extra attributes) 
#' where each row contains an edit restriction on purely categorical data. 
#' The function \code{editarray} converts (a vector of) edit(s)
#' in \code{character} or \code{expression} from to an \code{editarray} object. 
#' Edits may also be read from a \code{\link{data.frame}}, in which case it must have at least
#' a \code{character} column  with the name \code{edit}. It is not 
#' strictly necessary, but hightly recommended that the datamodel (i.e. the possible levels
#' for a variable) is included explicitly in the edits using an \code{\%in\%} statement, as shown
#' in the examples below. The function \code{\link{editfile}} can read categorical edits from a 
#' free-form text file.
#' 
#'
#'
#' @param editrules \code{character} or \code{expression} vector.
#' @param sep textual separator, to be used internally for separating variable from category names. 
#' @param env environment to evaluate the rhs of '==' or '\%in\%' in. 
#' @return \code{editarray} : An object of class \code{editarray}
#'
#' @example ../examples/editarray.R
#'
#' @seealso \code{\link{editrules.plotting}}, \code{\link{violatedEdits}}, \code{\link{localizeErrors}},
#'    \code{\link{editfile}}, \code{\link{editset}}, \code{\link{editmatrix}}, \code{\link{getVars}},
#'    \code{\link{blocks}}, \code{\link{eliminate}}, \code{\link{substValue}}, \code{\link{isFeasible}} 
#'    \code{\link{generateEdits}}, \code{\link{contains}}, \code{\link{is.editarray}}, \code{\link{isSubset}}
#'    
#' @export
editarray <- function(editrules, sep=":", env=parent.frame()){
    if (length(editrules) == 0 ) return(neweditarray(array(numeric(0),dim=c(0,0)),ind=list(),sep=sep))
    if ( is.data.frame(editrules) ){
        stopifnot("edit" %in% names(editrules))
        editrules <- editrules$edit
    }

    e <- parseEdits(editrules)
    v <- lapply(e,parseCat,sep=sep,env=env)
    
    # find always FALSE edits:
    iNull <- sapply(v,is.null)  
    
    # derive datamodel
    cols <- sort(unique(do.call(c,lapply(v[!iNull],names))))
    # get variable names
    vr <- sub(paste(sep,".*","",sep=""),"",cols)
    vars <- unique(vr)
    
    # get categories
    cat <- sub(paste(".+",sep,sep=""),"",cols)
    
    # build indexing list
    ind <- lapply(vars, function(v) which(v==vr))
    ind <- lapply(ind,function(I) {names(I) <- cat[I];I})
    names(ind) <- vars
    
    # edits with NA only extend the data model, is.null detects the always FALSE edit
    v <- v[sapply(v,function(u) is.null(u) || !is.na(u[1]))]
    # replace NULL with NA so the allways FALSE edit is included explicitly
    if ( length(v) > 0 )  v[sapply(v,is.null)] <- NA
        
    # set editarray values
    n <- length(cols)
    m <- length(v)
    if ( m == 0 ){
        return(neweditarray(
            array(logical(0),dim=c(m,n),dimnames=list(edits=NULL,variables=cols)),ind,sep)
        )

    }

    editnames <- paste("cat",1:m,sep="")
    E <- array(NA, dim=c(m,n), 
            dimnames = list(
                edits = editnames,
                variables = cols
            )
        )
    lapply(1:m,function(i) E[i,names(v[[i]])] <<- v[[i]])    
    # per variable, the boolean values not filled in during parsing must be derived.
    # they are the opposite from allready filled in values, or in case they are not involved,
    # all TRUE.
    for ( J in ind ){
        # vars not in any edit.
        I <- apply(E[,J,drop=FALSE],1,function(e) all(is.na(e)) ) 
        E[I,J] <- TRUE
        # vars in edits
        E[,J] <-  t(apply(E[,J,drop=FALSE],1,function(e){
            val <- e[!is.na(e)][1]
            e[is.na(e)] <- !val
            e
        }))
    }
    neweditarray(E,ind,sep=sep)
}

#' Derive textual representation from (partial) indices
#'
#' 
#'
#' @param ind a \code{list}, usually a (part of) the 'ind' attribute of an editarray
#' @param invert \code{logical} vector of  lenght length(ind)
#'
#' @return For every entry in \code{ind}, a character vector is returned where each entry
#'      is a statement of the form \code{"<var> \%in\% c('<cat1>',...,'<catN>')" or "<var> == '<catt>'"}
#'      if invert==TRUE. If invert==FALSE, the negation of the above statements is returned.
#'
#' 
#'
#' @keywords internal
ind2char <- function(ivd, ind=ivd, invert=logical(length(ivd)),useEqual=TRUE){
    v <- names(ivd)
    cats <- lapply(ivd, function(k) paste("'", names(k), "'", sep=""))
    op <- rep("%in%",length(ivd))
    l <- sapply(cats,length)
    if (useEqual){
        op[l == 1 & !invert] <- "=="
        op[l == 1 &  invert] <- "!="
    }
    cats[l>1] <- lapply(cats[l>1], function(cc) paste("c(",paste(cc,collapse=", "),")",sep=""))
    u <- paste(v,op,cats)
    u[l>1 &  invert] <- paste("!(",u[l>1 & invert],")")
    u <- sub("'FALSE'","FALSE",u)
    u <- sub("'TRUE'","TRUE",u)
    u <- sub("!= FALSE","== TRUE",u)
    u <- sub("!= TRUE","== FALSE",u)
    u
}


#' Convert to character
#'
#' @method as.character editarray
#' @param x editarray object
#' @param useIf \code{logical}. Use if( <condition> ) <statement> or !<condition> | <statement> ? 
#' @param datamodel \code{logical}. Include datamodel explicitly?
#' @param ... further arguments passed to or from other methods
#' @rdname editarray
#' @export
as.character.editarray <- function(x, useIf=TRUE, datamodel=TRUE, ...){
    A <- getArr(x)
    if (ncol(A) == 0 ){ 
        s <- character(nrow(A))
        if ( nrow(A)>0 ) names(s) <- rownames(A)
        return(s)
    }
    ind <- getInd(x)
    dm <- c()
    if ( datamodel ){
        dm <- ind2char(ind,useEqual=FALSE)
        names(dm) <- paste("dat",1:length(dm),sep="")
    }
    # edits
    if ( nrow(A) == 0 ) return(dm)
    edts <- character(nrow(A))
    for ( i in 1:nrow(A) ){
        a <- A[i,]
        involved <- sapply(ind, function(J) sum(a[J]) < length(J))
        # corner case: every record fails the edit (all categories TRUE).
        #if (!any(involved)) involved <- !involved
        ivd <- ind[involved]
        ivd <- lapply(ivd, function(J) J[a[J]])
        if ( length(ivd) == 0 ){
            edts[i] <- FALSE
        } else if ( length(ivd) == 1 ){
            edts[i] <- paste("if (",ind2char(ivd, ind),") FALSE")
        } else {
            
            n <- length(ivd)
            inv <- logical(n)
            inv[n] <- TRUE
            # corner case: all categories TRUE
            #if (all(involved) ) inv <- !inv
            ch <- ind2char(ivd, ind, invert=inv)
            if ( useIf ){
                edts[i] <- paste("if(", paste(ch[1:(n-1)],collapse=" & "), ")",ch[n])
# ALT:                edts[i] <- paste("if(", paste(ch,collapse=" & "), ") FALSE")
                
            } else {
                edts[i] <- paste("!(", paste(ch[1:(n-1)],collapse=" & "), ") |",ch[n])
# ALT:               edts[i] <- paste("!(", paste(ch,collapse=" & "), ") | FALSE")
            }
        }
    }
    names(edts) <- rownames(x)
    # add datamodel and return
    c(dm, edts) 
}


#' convert to data.frame
#' @method as.data.frame editarray
#' 
#' @rdname editarray
#' @return \code{as.data.frame}: \code{data.frame} with columns 'name', 'edit' and 'description'.
#' @export 
as.data.frame.editarray <- function(x, ...){
    edts <- as.character(x, ...)
    d <- data.frame(
        name=names(edts),
        edit=edts,
        row.names=NULL,
        stringsAsFactors=FALSE
    )
    if (!is.null(attr(x,'description'))) d$description <- attr(x,'description')
    d
}

#' Convert to expression 
#'
#' @rdname editarray
#' @export
#' @method as.expression editarray
as.expression.editarray <- function(x, ...){
  return(
    tryCatch(parse(text=as.character(x, ...)), 
             error=function(e){
               stop(paste("Not all edits can be parsed, parser returned", e$message,sep="\n"))
             }
             )
    )
}


#' editarray: logical array where every column corresponds to one
#' level of one variable. Every row is an edit. Every edit denotes
#' a *forbidden* combination.
#' @keywords internal
neweditarray <- function(E, ind, sep, names=NULL, levels=colnames(E),...){
    if ( is.null(names) & nrow(E)>0 ) names <- paste("cat",1:nrow(E),sep="")
    dimnames(E) <- list(edits=names,levels=levels)
    structure(E,
        class  = "editarray",
        ind    = ind,
        sep    = sep,
        ...
    )
}





