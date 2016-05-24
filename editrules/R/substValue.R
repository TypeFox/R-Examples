#' Replace a variable by a value in a set of edits.
#'
#' @note At the moment, objects of class \code{\link[=disjunct]{editenv}} are converted to \code{list}
#'  prior to processing (so no performance is gained there) and reconverted afterwards.
#'
#' @param E \code{\link{editset}}, \code{\link{editmatrix}}, \code{\link{editarray}}, 
#'      \code{\link[=disjunct]{editlist}} or \code{\link[=disjunct]{editenv}}
#' @param var \code{character} with name(s) of variable(s) to substitute
#' @param value vector with value(s) of variable(s)
#' @param ... arguments to be passed to or from other methods
#' @return \code{E}, with variables replaced by values
#' @example ../examples/substValue.R
#' @seealso \code{\link{eliminate}}
#' @export
#' @references
#'  Value substitution is extensively described in the package vignettes.
substValue <- function(E, var, value, ...){ 
    UseMethod("substValue")
}


# Reduce an editmatrix by substituting a variable
#
# Given a set of linear restrictions \eqn{E: {\bf Ax}\odot {\bf b}} with \eqn{\odot\in\{<,\leq,==\}},
# and matrix \eqn{{\bf A}} with columns \eqn{{\bf a}_1,{\bf a}_2,\ldots,{\bf a}_n}.
# Substituting variable \eqn{x_j} with a value \eqn{\tilde{\bf x}_j} means setting \eqn{{\bf a}_j=0}
# and \eqn{{\bf b}={\bf a}_j\tilde{x}_j}.
#
# Note that the resulting \code{\link{editmatrix}} may be inconsistent because of inconsistencies in
# \eqn{\tilde{\bf x}}.
#'
#' @method substValue editmatrix
#' @param reduce \code{logical} should the result be simplified? For \code{\link{editmatrix}} this has the same effect
#'  as calling the function \code{\link{reduce}}. For \code{\link{editarray}}, the datamodel of the substituted variable
#'  is reduced to a single value, and the variable itself is not removed. 
#' @param removeredundant \code{logical}. Should empty rows be removed?
#'
#' @rdname substValue 
#' @export
substValue.editmatrix <- function(E, var, value, reduce=FALSE, removeredundant=TRUE, ...){
  stopifnot(length(var)==length(value))
  if (length(var) == 0) return(E)
  
    v <- match(var, getVars(E), nomatch=0)
    if (any(v==0)){
        warning("Parameter var (", var[v==0], ") is not a variable of editmatrix E")
    }
    v <- v[v != 0]
    ib <- ncol(E)
   # typecast of 'value' so it may be passed as list (usefull in error localization).
    E[,ib] <- E[ ,ib] - E[ ,v]%*%as.numeric(value)

    if (reduce)
        E <- E[,-v, drop=FALSE]
    else 
        E[,v] <- 0
    if (removeredundant) {
        return( E[!isObviouslyRedundant.editmatrix(E),] )
    } else {
        return(E)
    }  
}



# Substitute a value in an editarray
#
# For an \code{\link{editarray}}, only rows with \code{<var>:<value>==TRUE} are kept.
# In the kept rows, categories not equal to <value> are set to \code{FALSE}
# If \code{reduce=TRUE}, columns corresponding to categories which are set
# to \code{FALSE} will be removed. Note that the function \code{\link{reduce}}
# has a different effect (it removes complete variables).
#
#' @method substValue editarray
#'
#'
#' @rdname substValue
#'
#' @export
substValue.editarray <- function(E, var, value, reduce=FALSE, ...){
  stopifnot(length(var)==length(value))
  if (length(var) == 0) return(E)
  
  ind <- getInd(E)
    sep=getSep(E)
    A <- getArr(E)
    value <- as.character(value)
    for ( i in seq_along(var) ){
        vr <- var[i]
        vl <- value[i]
        J <- ind[[vr]]
        ii <- J[vl]
        if ( is.null(ii) || is.na(ii) ) 
            stop(paste("Variable ", vr,"not present in editarray or cannot take value",vl))

        I <- A[,ii]
        if ( reduce ){
            A <- A[ ,-setdiff(J,ii) ,drop=FALSE]
            ind <- indFromArray(A, sep)
        } else {
            A[,J] <- TRUE
        }
    }
    neweditarray(
        E = A[I,,drop=FALSE], 
        ind = ind, 
        sep = sep, 
        levels = colnames(A) 
    )
}

#' Compute index from array part of editarray
#' 
#' @param A boolean array
#' @param sep separator
#' @keywords internal
#'
indFromArray <- function(A,sep){
    if (ncol(A) == 0 ) return(list())
    cn <- colnames(A)
    l <- strsplit(cn,sep)
    V <- sapply(l,`[`,1)
#    C <- sapply(l,`[`,-1)
    C <- sapply(l,function(g) ifelse(length(g[-1])==1,g[-1],""))
    vars <- unique(V)
    ind <- lapply(vars, function(v) which(v==V))
    names(ind) <- vars
    ind <- lapply(ind, function(k){ names(k) <- C[k]; k})
    ind
}



# Substitute values in an \code{\link{editset}}
#
# For an \code{\link{editset}}, purely numerical variables are
# substitutes as in an \code{\link{editmatrix}} and categorical
# as in an \code{\link{editarray}}. Numerical variables appearing
# logical constraints are substituted and if truth values can
# be derived these are substituted in the logical constraint.
# 
#' @param simplify Simplify editset by moving logical edits containing a single
#'      numerical statement to the pure numerical part? (This is mostly for internal purposes
#'      and overwriting the default should normally not be necessary for package users).
#'
#' @method substValue editset
#' 
#' @rdname substValue 
#' @export
substValue.editset <- function(E, var, value, simplify=TRUE, ...){
# Techical note. Substituting a dummy variable (e.g. .l1) with TRUE or 
# FALSE amounts to making an assumption about the validity
# of the condition stated in that dummy. As such, it should not be added
# to the numerical editmatrix (since that editmatrix is only relevant when the
# assumed condition is already fulfilled). Instead, the condition is 
# added to the 'condition' attribute of the editset.
  
#TODO make it possible to supply value = list(x=1, A="a") which makes
# substituting values in an editset a lot easier. Especially when we have 
# used localizeErrors and want the solution space by substituting non adapted variables.
  stopifnot(length(var)==length(value))
  if (length(var) == 0) return(E)
  
    catidx <- var %in% getVars(E, type="cat")

    # the nonnumeric case is simple
    if ( !is.numeric(value) ){
        E$mixcat <- substValue(E$mixcat,var,value,...)
        # move substituted dummies to "condition"
        id <- var %in% getVars(E,type='dummy')
        if ( any(id) ){
            dvar <- rownames(E$mixnum) %in% var[id]
            v <- as.character(E$mixnum[dvar,])
            v[!value[id]] <- invert(v[!value[id]])
            attr(E,"condition") <- c(editmatrix(v),attr(E,"condition"))
            E$mixnum <- E$mixnum[!dvar,]
        }
        if ( simplify ) E <- simplify(E)
        return(E)
    }
    # substitute pure numeric data
    i1 <- var %in% getVars(E$num)
    if ( any(i1) ){ # time-saving condition
        numvar <- var[i1]
        numval <- value[i1]
        innum <- colSums(contains(E$num, numvar )) > 0
        if ( any(innum) ) 
            E$num <- substValue(E$num, numvar[innum], numval[innum])
    }
    # substitute in condition 
    cnd <- condition(E)
    if ( var %in% getVars(cnd) ) condition(E) <- substValue(cnd,var,value,...)
    # substitute in then-clauses
    i1 <- var %in% getVars(E$mixnum)
    if ( any (i1) ){ # time-saving condition
        mixvar <- var[i1]
        mixval <- value[i1]
        u <- contains(E$mixnum, mixvar)
        inmix <- colSums(u) > 0
        if ( any(inmix) ){
            E$mixnum <- substValue(
                E$mixnum, 
                mixvar[inmix], 
                mixval[inmix],
                removeredundant=FALSE
            )
            # did substitution yield any certainties?
            cntr <- isContradiction(E$mixnum)
            taut <- isTautology(E$mixnum)
            # dummy variables to be eliminated from mixcat
            lvar <- apply(u[,inmix,drop=FALSE],1,any)
            dvar <- rownames(u)
            dval <- logical(length(taut))
            dval[lvar & cntr] <- FALSE
            dval[lvar & taut] <- TRUE
            isub <- lvar & (cntr | taut)
            E$mixcat <- substValue(E$mixcat, dvar[isub],dval[isub])
        }
    }
    if ( simplify ) E <- simplify(E)
    removeRedundantDummies(E)
}




# Returns which linear edits are obvious contradictions.
#  - Accurate to 8 figures.
#  - Assumes editmatrix normality
isContradiction <- function(E){
    tol = 1e-8
    ops <- getOps(E)
    absA <- abs(getA(E))
    nil <- rowSums(absA) < ncol(absA)*tol
    b <- getb(E)
    I <- logical(nrow(absA))
    eq <- ops=='=='
    lt <- ops=='<'
    le <- !eq & !lt
    I[eq] <- nil[eq] & abs(b[eq]) > tol
    I[lt] <- nil[lt] & b[lt] <= 0 
    I[le] <- nil[le] & b[le] < tol
    I
}

# returns which linear edits are obviously TRUE
# - Accurate to 8 figures
# - Assumes editmatrix normality
isTautology <- function(E, tol=sqrt(.Machine$double.eps)){
    tol = 1e-8
    ops <- getOps(E)
    absA <- abs(getA(E))
    nil <- rowSums(absA) < ncol(absA)*tol
    b <- getb(E)
    I <- logical(nrow(absA))
    eq <- ops=='=='
    lt <- ops=='<'
    le <- !eq & !lt
    I[eq] <- nil[eq] & abs(b[eq]) < tol
    I[lt] <- nil[lt] & b[lt] > tol
    I[le] <- nil[le] & b[le] >= -tol
    I
}



#'
#' @method substValue editlist
#' @rdname substValue
#' @export
substValue.editlist <- function(E, var, value, ...){
    L <- varTypeAndOccurrence(E,var)
    if ( length(L) == 1 && is.na(L) ){
        return(E)      
    }
    type = L$type
    iRemove <- logical(length(E))
    for ( i in which(L$occurs) ){
        E[[i]] <- substValue(E[[i]],var,value,...)
        if ( !isFeasible(condition(E[[i]])) ) iRemove[i] <- TRUE
    }
    E[!iRemove]
}


#' @method substValue editenv
#' @rdname substValue
#' @export
substValue.editenv <- function(E,var,value,...){
    L <- as.list(E)
    L <- substValue.editlist(E)
    list2env(L)
}



