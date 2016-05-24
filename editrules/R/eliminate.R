
#' Eliminate a variable from a set of edit rules
#' 
#' Eliminating a variable amounts to deriving all (non-redundant) edits not containing
#' that variable. Geometrically, it can be seen as a projection of the solution space
#' (records obeying all edits) along the eliminated variable's axis. If the solution space
#' is non-concex (as is the usually case when conditional edits are involved), multiple
#' projections of convex subregions are performed.
#'
#' @param E \code{\link{editmatrix}} or \code{\link{editarray}} 
#' @param var name of variable to be eliminated
#' @param ... argumemts to be passed to or from other methods
#' @export
#' @seealso \code{\link{substValue}}, \code{\link{isObviouslyInfeasible}}, \code{\link{isObviouslyRedundant}},
#'  \code{\link{generateEdits}}
#'
#' @return If \code{E} is an \code{\link{editmatrix}} or \code{\link{editarray}}, an object of the same class
#'   is returned. A returned \code{\link{editmatrix}} contains an extra \code{history} attribute which is used
#'   to reduce the number of generated edits in consecutive eliminations (see \code{\link{getH}}). If \code{E} is an \code{\link{editset}},
#'   an object of class \code{\link[=disjunct]{editlist}} is returned. 
eliminate <- function(E, var, ...){
    UseMethod("eliminate")
}


#' For objects of class \code{\link{editmatrix}}, Fourier-Motzkin elimination 
#' is used to eliminate a variable from the of linear (in)equality restrictions.
#' An observation of Kohler (1967) is used to reduce the number of implied 
#' restrictions. Obvious redundancies of the type 0 < 1 are removed as well. 
#'
#' @method eliminate editmatrix
#' 
#'
#'
#' @rdname eliminate 
#' @example ../examples/eliminate.R
#'
#' @references
#' D.A. Kohler (1967) Projections of convex polyhedral sets, Operational Research
#' Center Report , ORC 67-29, University of California, Berkely.
#' 
#' H.P. Williams (1986) Fourier's method of linear programming and its dual,
#' The American Mathematical Monthly 93, 681-695
#' @export
eliminate.editmatrix <- function(E, var, ...){
    if (!isNormalized(E)) E <- normalize(E)
    vars <- getVars(E)

    if (!var %in% vars){
        warning("Trying to eliminate variable not occuring in E")
        return(E)
    }
    d <- getH(E)   
    n <- geth(E)  
    if ( is.null(d) ){
        n <- 0
        d <- matrix(FALSE,nrow=nrow(E),ncol=nrow(E))
        diag(d) <- TRUE 
        colnames(d) <- rownames(E)
    }
    
    
    m <- as.matrix(E)
    ops <- getOps(E)
    
    coefs <- m[,var]
    I <- coefs != 0
    
    eq <- I & (ops == "==")
    
    upper <- which(!eq & coefs > 0)
    lower <- which(!eq & coefs < 0)
    
    coefs[!eq] <- abs(coefs[!eq])
    
    eq <- which(eq);

    # elimination possible?
    if ( (length(upper) > 0 && length(lower) > 0) ||
         (length(eq) >= 1 && (length(upper) > 0 || length(lower) > 0)) ||
         (length(eq) >= 2) ){
       n <- n+1
    } else {
        return(E[ E[ ,var]==0,])
    } 


    #normalize matrix, every row has coefficient 1 or -1 for var
    m[I,] <- m[I,] / coefs[I]

    # eqs and ineqs w/coef>0 ==> ineqs w/coef<0
    equpper <- c(eq, upper)
    I1 <- rep(equpper,each=length(lower))
    I2 <- rep(lower, times=length(equpper))
    ml <- m[I1,,drop=FALSE] + m[I2,,drop=FALSE]
    ol <- ifelse(ops[I1] != "<", ops[I2], ops[I1])
    dl <- d[I1,,drop=FALSE] | d[I2,,drop=FALSE]

    # eqs ==> ineqs w/coef>0
    I1 <- rep(eq,each=length(upper))
    I2 <- rep(upper, times=length(eq))
    mu <- m[I2,,drop=FALSE] - m[I1,,drop=FALSE]
    ou <- ops[I2]
    du <- d[I1,,drop=FALSE] | d[I2,,drop=FALSE]

    # eqs ==> eqs
    me <- m[logical(0),,drop=FALSE]
    de <- d[logical(0),,drop=FALSE]
    if ( length(eq)>1){
        me <- t(t(m[eq[-1],,drop=FALSE]) - m[eq[1],])
        de <- t(t(d[eq[-1],,drop=FALSE]) | d[eq[1],])       
    } 
    oe <- rep("==",nrow(me))

    m <- rbind(ml,mu,me,m[!I,,drop=FALSE])
    d <- rbind(dl,du,de,d[!I,,drop=FALSE])
    o <- c(ol,ou,oe,ops[!I])
    redundant <- rowSums(d) > n + 1 | isObviouslyRedundant.matrix(E=m, operators=o)

    m <- m[!redundant,,drop=FALSE]
    d <- d[!redundant,,drop=FALSE]
    if ( nrow(m) > 0 ){
      rownames(m) <- paste("num",1:nrow(m),sep="")
    }
    neweditmatrix(
          m
        , o[!redundant]
        , normalized = TRUE 
        , H = d 
        , h = n
    )
}




#' For categorical edits in an \code{\link{editarray}}, the elimination 
#' method is based on repeated logical reduction on categories. See Van der Loo (2012) for 
#' a description. 
#'
#' @method eliminate editarray
#' @rdname eliminate
#' @references
#'
#' M.P.J. van der Loo (2012) Variable elimination and edit generation with a flavour of 
#'  semigroup algebra (submitted).
#'
#' @export
eliminate.editarray <- function(E, var, ...){
    # do not bother with edits not containing var
    I <- contains(E,var)
    # nothing to eliminate...
    if ( sum(I) == 0 ) return(E) 

    ind <- getInd(E)
    J <- ind[[var]]

    # if elimination is not possible... (at least one category cannot be resolved)
    if ( any(colSums(E[I,J,drop=FALSE])==0) ) return(E[!I,])

    A <- getArr(E)
    At <- A[!I,,drop=FALSE]
    A <- A[I,,drop=FALSE]


    k <- integer(0)
    for ( j in J ){
        if (nrow(A) <= 1) break
        aPlus <- A[ A[,j],,drop=FALSE]
        aMin  <- A[!A[,j],,drop=FALSE]  
        m <- nrow(aMin)
        n <- nrow(aPlus)
        if ( m == 0 ) next
        if ( n == 0 ) break
        B <- array(FALSE,dim=c(n*m,ncol(A)))
        I1 <- rep(1:n,times=m)
        I2 <- rep(1:m,each=n)
        B[I1,J] <-  aPlus[I1,J,drop=FALSE] | aMin[I2,J,drop=FALSE]
        B[I1,-J] <- aPlus[I1,-J,drop=FALSE] & aMin[I2,-J,drop=FALSE]
        B <- rbind(B,aPlus)
        A <- B[!isRedundant.boolmat(B,ind),,drop=FALSE]
        A <- A[!isSubset.boolmat(A),,drop=FALSE]
        el <- apply(A[,J,drop=FALSE],1,all)
        At <- rbind(At,A[el,,drop=FALSE])     
        A <- A[!el,,drop=FALSE]
    }

    neweditarray(At,ind=ind,sep=getSep(E),levels=getlevels(E))
}

#'
#' For an \code{\link{editset}}, \code{E} is transformed to an \code{\link[=disjunct]{editlist}}.
#' Each element of an \code{\link[=disjunct]{editlist}} describes a convex subregion of the 
#' total solution space of the \code{\link{editset}}. After this, the elimination method for 
#' \code{\link[=disjunct]{editlist}} is called.
#'
#' @export 
#' @rdname eliminate
#' @method eliminate editset
eliminate.editset <- function(E,var,...){
    D <- disjunct(E)
    eliminate.editlist(D,var,...)
}


#'
#' For an \code{\link[=disjunct]{editlist}}, the variable is eliminated from each
#' consituting \code{\link{editset}}. 
#'
#' @export
#' @rdname eliminate
#' @method eliminate editlist
#'
eliminate.editlist <- function(E, var, ...){
    # determine variable type (num or cat)
    L <- varTypeAndOccurrence(E,var)
    if ( length(L) == 1 && is.na(L) ){
        return(E)
    } else {
        type <- L$type
        for ( i in which(L$occurs) ){
            E[[i]][[type]] <- eliminate(E[[i]][[type]],var)
        }
    }
    E
}




