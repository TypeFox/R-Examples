#' Read general edits
#' 
#' An \code{editset} combines numerical (linear), categorical and conditional restrictions
#' in a single object. Internally, it consists of two \code{\link[=editmatrix]{editmatrices}} 
#' and an \code{\link{editarray}}.
#'
#' The function \code{editset} converts a \code{character} or \code{expression} vector to an editset.
#' Alternatively, a \code{data.frame} with a column called \code{edit} can be supplied. Function
#' \code{\link{editfile}} reads edits from a free-form textfile.
#'
#' 
#'
#' @param editrules \code{character} vector, \code{expression} vector or \code{data.frame} (see details) containing edits.
#' @param env environment to parse categorical edits in (normally, users need not specify this)
#'
#' @seealso \code{\link{editrules.plotting}}, \code{\link{violatedEdits}}, \code{\link{localizeErrors}},
#'  \code{\link{getVars}}, \code{\link{disjunct}}, \code{\link{eliminate}}, \code{\link{substValue}}, 
#'  \code{\link{isFeasible}}, \code{\link{contains}}, \code{\link{is.editset}}
#' @example ../examples/editset.R
#' @return \code{editset}: An object of class \code{editset}
#' @export
editset <- function(editrules, env=new.env()){
    if (is.data.frame(editrules)){
        stopifnot("edit" %in% names(editrules))
        editrules <- editrules$edit
    }
    
    # detect edit types
    num <- parseEdits(editrules, type="num")
    cat <- parseEdits(editrules, type="cat")
    mix <- parseEdits(editrules, type="mix")
  

    # pure numerical edits    
    num <- editmatrix(num)
    nnum <- nrow(num)
    # pure categorical edits    
    cat <- editarray(cat,env=env)
    ncat <- nrow(cat)
    if ( ncat > 0 ) rownames(cat) <- paste("cat",(nnum+1):(nnum+ncat),sep="")


    nmix <- length(mix)
    mixnum <- expression()
    mixcat <- vector(mode="expression", nmix)
    if ( nmix > 0 ){ 
        mixl <- vector(mode="list", nmix)
        nms <- names(mix)
        numid <- 0
  
        for (i in seq_along(mix)){
            m <- parseMix(mix[[i]], nms[i], numid=numid)
            numid <- m$numid
            mixl[[i]] <- m
            mixcat[[i]] <- m$cat
            mixnum <- c(mixnum, m$nums)
        }

        # combine categorical edits, datamodel with dummies for mixed edits and mixed edits
        mix <- c(
            ind2char(getInd(cat)),
            if ( length(mixnum) > 0 ) paste(names(mixnum), "%in% c(TRUE,FALSE)"),
            as.character(mixcat)
        )
    }

    # stick categorical and mixed edits in one editarray
    mixcat <- editarray(mix, env=env)
    nmix <- nrow(mixcat)
    if ( nmix>0 ){
        rownames(mixcat) <- paste("mix",(nnum+ncat+1):(nnum+ncat+nmix),sep="")
        mixcat <- c(cat, mixcat)
    } else {
        mixcat <- cat
    }

    # create editmatrix for mixed edits and name them with dummy variable names
    nms <- names(mixnum)
    mixnum <- editmatrix(as.character(mixnum))
    rownames(mixnum) <- nms

    removeRedundantDummies( 
        neweditset(
            num=num,
            mixcat=mixcat,
            mixnum=mixnum
        )
    )
}

#
#
#
neweditset <- function(num, mixcat, mixnum, condition=editmatrix(expression()), ...){

  structure(
      list( num = num
          , mixnum = mixnum
          , mixcat = mixcat
          ),
    class="editset", 
    condition=condition,
    ...
  )

}

#' Get condition matrix from an editset. 
#'
#' @param E an \code{\link{editset}}
#' @return an \code{\link{editmatrix}}, holding conditions under which the editset is relevant.
#' @export
#' @seealso \code{\link{disjunct}}, \code{\link{separate}}, \code{\link{editset}}
condition <- function(E) attr(E,'condition')

`condition<-` <- function(x, value){
    if (!is.editset(x) ) stop("only for editset")
    attr(x,'condition') <- value
    x
}



#' Add dummy variable to the data.frames, these are needed for errorlocations etc.
#' @param E editset
#' @param dat data.frame
#' @return data.frame with dummy variables added
#' @keywords internal
adddummies <- function(E, dat){
  dummies <- !violatedEdits(E$mixnum, dat)
  cbind(dat, dummies)
}



#' @method as.character editset
#'
#' @param x an \code{\link{editset}}
#' @param datamodel include datamodel?
#' @param useIf return vectorized version?
#' @param dummies return datamodel for dummy variables?
#' @param ... arguments to be passed to or from other methods
#' @rdname editset
#' @export
as.character.editset <- function(x, datamodel=TRUE, useIf=TRUE, dummies=FALSE, ...){
    num <-  as.character(x$num)
    numc <- as.character(x$mixnum)
    catc <- as.character(x$mixcat, datamodel=datamodel, useIf=useIf)
    for ( n in names(numc) ){
        catc <- catc[!grepl(paste(n,'%in%'),catc)]
        catc <- gsub(paste(n,'== FALSE'), invert(numc[n]),catc)
        catc <- gsub(paste(n,'== TRUE'), numc[n],catc)
    }
    # remove datamodel which are a consequence of conditional edits
    
    c(num,catc)
}


# invert a textual numerical edit
invert <- function(e){
    gte   <- grepl(">=",e)
    gt    <- grepl(">",e) & !gte
    lte   <- grepl("<=",e)
    lt    <- grepl("<",e) & !lte
    eq    <- grepl("==",e)
    ineq  <- grepl("!=",e)
    e[gte]  <- gsub(">=","<",e[gte])
    e[gt]   <- gsub(">","<=",e[gt])
    e[lte]  <- gsub("<=",">",e[lte])
    e[lt]   <- gsub("<",">=",e[lt])
    e[ineq] <- gsub("!=","==",e[ineq])
    e[eq]   <- gsub("==","!=",e[eq])
    e
}



#'
#' @method as.data.frame editset
#' @rdname editset
#' @export 
#' @return \code{as.data.frame}: a \code{data.frame} with columns 'name' and 'edit'. 
#   If the editset has a \code{description} attribute, a third column named 'description' is added.
as.data.frame.editset <- function(x, ...){
    edts <- as.character(x, datamodel=TRUE,...)
    d <- data.frame(
        name=names(edts),
        edit=edts,
        row.names=NULL,
        stringsAsFactors=FALSE
    )
    if (!is.null(attr(x,'description'))) d$description <- attr(x,'description')
    d
}




#' Determine edittypes in editset based on 'contains(E)'
#'
#' Determines edittypes based on the variables they contain (not on names of edits).
#'
#' @param E editset
#' @param m if you happen to have \code{contains(E)} handy, it needs not be recalculated.
#' @seealso \code{\link{contains}}
#' @export
editType <- function(E, m=NULL){
    if ( !is.editset(E) ) stop('Argument is not of class editset')
    type <- vector(length=nedits(E),mode='character')
    nnum <- nrow(E$num)
    if ( is.null(m) ) m <- contains(E)
    if ( nnum > 0 ){
        type[1:nnum] <- 'num'
        m <- m[-(1:nnum),,drop=FALSE]
    }

    catvar <- getVars(E,'cat')
    mixvar <- getVars(E,'mix')

    icat <- rowSums(m[,catvar,drop=FALSE])>0
    imix <- rowSums(m[,mixvar,drop=FALSE])>0
    type[which(imix) + nnum] <- 'mix'
    type[which(icat&!imix) + nnum] <- 'cat'
    type
}


#' Coerce x to an editset
#'
#' \code{x} may be an editset, editmatrix, editarray or character vector
#' @param x object or vector to be coerced to an editset
#' @param ... extra parameters that will be passed to \code{as.character}, if necessary
#' @export
as.editset <- function(x, ...){
  if (is.editset(x))
    return(x)
  
  if (is.cateditmatrix(x)){
    return (editset(as.character(x, asIfStatement=TRUE)))
  }
  
  if (is.editmatrix(x)){
    E <- editset(expression())
    E$num <- x
    return(E)
  }
  
  if (is.editarray(x)){
    E <- editset(expression())
    E$mixcat <- x
  }
  
  editset(as.character(x, ...))
}


#' Simplify logical mixed edits in an editset
#'
#' Logical edits consisting of a single numerical statement are
#' inverted and added to the \code{editmatrix} (\code{\$num}) of the editset.
#'
#' @param E an editset
#' @param m \code{contains(E)}. Speeds up calculation if you have one handy.
#'
#' @keywords internal
simplify <- function(E, m=NULL){
    # NOTE: might be interesting for @export
    if (!is.editset(E)) stop("Argument not of class 'editset'")
    mixvar <- getVars(E,type='mix')
    catvar <- getVars(E,type='cat')
    dummies <- getVars(E,type='dummy')
    if ( is.null(m) ) m <- contains(E) 
    
    # move pure numerical edits to $num
    g <- contains(E$mixcat)
    r <- rowSums(g[,dummies,drop=FALSE]) == 1 & rowSums(g[,catvar,drop=FALSE]) == 0

    if ( any(r) ){ 
        # this is damned ugly, but: convert to character
        v <- as.character(E[which(r)+nrow(E$num),,drop=FALSE],datamodel=FALSE)
        # strip if ( ) FALSE to make parsing to editmatrix possible
        v <- sub("\\).+$", '',sub("^.+\\(",'',v))  
        v <- invert(v)
        E$mixcat <- reduce(E$mixcat[!r,,drop=FALSE])
        E$num <- c(E$num,editmatrix(v))
        E$mixnum <- reduce(E$mixnum[dummies[dummies %in% getVars(E$mixcat)],])
    } else {
        E$mixcat <- reduce(E$mixcat)
        E$mixnum <- reduce(E$mixnum[dummies[dummies %in% getVars(E$mixcat)],])
    }
   
    E
}   


#' Remove redundant dummy variables
#'
#' Remove duplicated dummy variables from an editset.
#'
#' @param E \code{\link{editset}}
#' @param tol positive number
#' @keywords internal
#'
removeRedundantDummies <- function(E, tol=1e-8){
    if ( is.null(E$mixnum) || nrow(E$mixnum) < 2) return(E)
    if (!isNormalized(E$mixnum)) E$mixnum <- normalize(E$mixnum)
    
    op2num <- c("<" = 1, "<=" = 2, "==" = 3, ">=" = 4, ">" = 5)

    d <- dist(cbind(getAb(E$mixnum), op2num[getOps(E$mixnum)]))
    id <- d < tol
    if ( !any( id ) ) return(E)
    m <- nrow(E$mixnum)
    dupnames <- rownames(E$mixnum)[2:m]
    orgnames <- rownames(E$mixnum)[1:(m-1)]
    v <- matrix(FALSE,nrow=m-1, ncol=m-1, dimnames=list(dupnames,orgnames))
    v[lower.tri(v,1)] <- d < tol
    v <- v[apply(v,1,any),,drop=FALSE]
    v <- v[,!apply(!v,2,all),drop=FALSE]
    v <- v[,!colnames(v) %in% rownames(v),drop=FALSE]

    # remove duplicates from mixnum
    dupvars <- rownames(v)
    w <- rownames(E$mixnum)
    E$mixnum <- E$mixnum[!w %in% dupvars,]

    # replace original values by duplicates in mixcat; remove duplicates
    A <- getArr(E$mixcat)
    ind <- getInd(E$mixcat)
    sep <- getSep(E$mixcat)

    w <- which(v,arr.ind=TRUE)
    dups <- rownames(v)[w[,1]]
    names(dups) <- colnames(v)[w[,2]]
    while( length(dups) > 0 ){
        dup <- dups[1]
        org <- names(dup)
        iM <- contains.boolmat(A, ind, dup)
        A[iM,ind[[org]]] <- A[iM,ind[[dup]],drop=FALSE]
        dups <- dups[-1]
    }

    idup <- unlist(ind[dupvars])
    A <- A[,-idup,drop=FALSE]
    ind <- indFromArray(A,sep)
    E$mixcat <- neweditarray(A,ind=ind,names=rownames(A),sep=sep)

    E
}



## quick test
# es <- editset(expression(if (x > 0) y + 1 < 20
#                         , x <= 100
#                         , if (x < 10) y >= 2
#                         , A %in% c("a1", "a2")
#                         , B %in% c("b1", "b2")
#                         , if (A == "a1") B == "b2"
#                         , if (y > 0) A == "a2"
#                         )
#              )
# #es
# 
# dat <- data.frame(x=1:2, y=10:9, A="a1", B="b2")
# #adddummies(es,dat)
# 
# violatedEdits(es, dat)
