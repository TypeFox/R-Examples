# This is file ../spam/R/spam_solve.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     

########################################################################
########################################################################
#
# Contains routines linked to solving spd linear systems. Namely:
#    chol, solve, backsolve, forwardsolve
#    determinant                (the later because it is based on chol)
#
# As well as associated S4 elements.
# 
# The key element is a new class: "spam.chol.NgPeyton", the output
# of 'chol'
# 
########################################################################
########################################################################

setClass("spam.chol.NgPeyton",
         representation(entries="numeric",      colindices="integer",
                        colpointers="integer",  rowpointers="integer",
                        dimension="integer",
                        pivot="integer",        invpivot="integer",
                        supernodes="integer",   snmember="integer",
                        memory="integer",       nnzA="integer")
         )

# lindx=  colindices
# xlindx= colpointers
# xlnz=   rowpointers
# snode=snmember
# xsuper=supernodes
# c(... nnztmp,cachesize)= memory

  
#setClass("spam.chol.NgPeyton",
#         representation(nrow="integer",nnzlindx="integer",
#                        nsuper="integer",lindx="integer",xlindx="integer",nnzl="integer",
#                        lnz="numeric",xlnz="integer",invp="integer",perm="integer",
#                        xsuper="integer"),
         # the prototype corresponds to the cholesky of '1'
#         prototype=prototype(nrow=as.integer(1),nnzlindx=as.integer(1),
#           nsuper=as.integer(1),lindx=as.integer(1),xlindx=as.integer(c(1,2)),
#           nnzl=as.integer(1),lnz=1.0,xlnz=as.integer(c(1,2)),
#           invp=as.integer(1),perm=as.integer(1),xsuper=as.integer(c(1,2))
#           )
#         )


########################################################################

print.spam.chol.NgPeyton <- function(x,...) {
  nrow <- x@dimension[1]
  nnzR <- x@rowpointers[nrow+1]-1
  cat("(Upper) Cholesky factor of dimension ", nrow,
                "x", nrow, " with ",nnzR," nonzero elements.", sep = "", fill=TRUE)
  cat("    (The object is supposed to be used with: 'as.spam', 'backsolve', 'forwardsolve', etc.)\n",
      fill=TRUE)
  cat("Class 'spam.chol.NgPeyton'\n")
  invisible(NULL)
}
summary.spam.chol.NgPeyton <- function(object,...) {
  nrow <- object@dimension[1]
  nnzR <- object@rowpointers[nrow+1]-1
  dens <- nnzR/(nrow^2)
  nnzc <- length(object@colindices)
  fill <- nnzR/((object@nnzA+nrow)/2)
  cat("(Upper) Cholesky factor of class 'spam.chol.NgPeyton' of dimension ", nrow,
                "x", nrow, " with ",nnzR," (row-wise) nonzero elements.", sep = "", fill=TRUE)
  cat("    Density of the factor is ", signif(dens * 100, 3),"%.\n", sep = "")
  cat("    Fill-in ratio is ", signif(fill, 3),"\n", sep = "")
  cat("    (Optimal argument for 'chol' is 'memory=list(nnzR=",nnzR,
            ifelse(object@nnzA<nnzc,paste(",nnzcolindices=",nnzc, sep = ""),""),")'.)\n", sep = "")
  cat("Class 'spam.chol.NgPeyton'\n")
  invisible(list(nnzR=nnzR,nnzcolindices=nnzc,density=dens,fillin=fill))
}

setMethod("show","spam.chol.NgPeyton", function(object) {
  nrow <- object@dimension[1]
  nnzR <- object@rowpointers[nrow+1]-1
  cat("(Upper) Cholesky factor of dimension ", nrow,
                "x", nrow, " with ",nnzR," (row-wise) nonzero elements.", sep = "", fill=TRUE)
  cat("    (The object is supposed to be used with: 'as.spam', 'backsolve', 'forwardsolve', etc.)\n",
      fill=TRUE)
  cat("Class 'spam.chol.NgPeyton'\n")
  invisible(NULL)
        })

"diag.of.spam.chol.NgPeyton" <- function(x, nrow, ncol)
  return( x@entries[x@rowpointers[-(x@dimension[1]+1)]])


setMethod("diag",    "spam.chol.NgPeyton", diag.of.spam.chol.NgPeyton)
#setMethod("diag<-",  "spam.chol.NgPeyton", function(x,...) stop("operation not allowed on 'spam.chol.NgPeyton' object"))
setMethod("print",   "spam.chol.NgPeyton", print.spam.chol.NgPeyton)
setMethod("summary", "spam.chol.NgPeyton", summary.spam.chol.NgPeyton)
setMethod("dim",     "spam.chol.NgPeyton",function(x) x@dimension)
setMethod("length",  "spam.chol.NgPeyton",function(x) x@rowpointers[x@dimension[1]+1]-1)
setMethod("length<-","spam.chol.NgPeyton",function(x,value) stop("operation not allowed on 'spam.chol.NgPeyton' object") )
setMethod("dim<-",   "spam.chol.NgPeyton",function(x,value) stop("operation not allowed on 'spam.chol.NgPeyton' object") )

setMethod("c","spam.chol.NgPeyton", function(x,...,recursive=TRUE){
  nrow <- x@dimension[1]
  nnzR <- x@rowpointers[nrow+1]-1
  newx <- new("spam")
  nsuper <- as.integer( length(x@supernodes)-1)
  xcolindices <- .Fortran('calcja',
                           nrow, nsuper, x@supernodes, x@colindices, x@colpointers, x@rowpointers,
                           xja=vector("integer",nnzR),
                           NAOK = .Spam$NAOK,PACKAGE = "spam")$xja
  cx <- .Fortran("spamcsrdns",
                 nrow=nrow,
                 entries=as.double(x@entries),
                 colindices=xcolindices,
                 rowpointers=x@rowpointers,
                 res=vector("double",nrow*nrow),  
                 NAOK=.Spam$NAOK,PACKAGE = "spam")$res
  if (length( list(...)) < 1)
    return( cx)
  else
    c( cx,c(...,recursive),recursive)
})


"as.spam.chol.NgPeyton" <- function(x, eps = .Spam$eps){
  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
  nrow <- x@dimension[1]
  nnzR <- x@rowpointers[nrow+1]-1
  newx <- new("spam")
  nsuper <- as.integer( length(x@supernodes)-1)
  slot(newx,"entries",check=FALSE) <- x@entries
  slot(newx,"colindices",check=FALSE) <- .Fortran('calcja',
                           nrow, nsuper, x@supernodes, x@colindices, x@colpointers, x@rowpointers,
                           xja=vector("integer",nnzR),
                           NAOK = .Spam$NAOK,PACKAGE = "spam")$xja
  slot(newx,"rowpointers",check=FALSE) <- x@rowpointers
  slot(newx,"dimension",check=FALSE) <- x@dimension
  return(newx)
}

setMethod("as.spam","spam.chol.NgPeyton", as.spam.chol.NgPeyton)


##NB setGeneric("backsolve", def = function(r, x, ...) standardGeneric("backsolve"),
#                 useAsDefault= function(r, x,...) base::backsolve(r, x, ...))

# We have some issues here... hence I postpone the proper implementation!!!
# http://r.789695.n4.nabble.com/class-extension-and-documentation-tt4161373.html#none

#"backsolve" <- function(r,x, ...) UseMethod("backsolve")
#"backsolve.default" <- base::backsolve
#setGeneric("backsolve")
#setMethod("backsolve","matrix",base::backsolve)
  
#"forwardsolve" <- function(l,x, ...) UseMethod("forwardsolve")
#"forwardsolve.default" <- base::forwardsolve
#setGeneric("forwardsolve")
#setMethod("forwardsolve","matrix",base::forwardsolve)


setGeneric("backsolve", def = function(r, x, ...) standardGeneric("backsolve"),
           useAsDefault= function(r, x, ...) base::backsolve(r, x, ...))

setGeneric("forwardsolve", def = function(l, x, ...) standardGeneric("forwardsolve"),
           useAsDefault= function(l, x, ...) base::forwardsolve(l, x, ...))

# adapted from methods
#setGeneric("forwardsolve", function(l, x, k, upper.tri = FALSE, transpose = FALSE, ...)
#           standardGeneric("forwardsolve"),
#           useAsDefault = function(l, x, k = ncol(l), upper.tri = FALSE, transpose = FALSE, ...)
#                  base::forwardsolve(l, x, k = k, upper.tri = upper.tri, transpose = transpose, ... ),
#           signature = c("l", "x"))#, where = where)
##### setGenericImplicit("forwardsolve")#, restore=FALSE)

  
"ordering.default" <- function(x,inv=FALSE) stop('Operation not defined form this class')

#ordering <- function(x,...) stop('Operation not defined form this class')
#setGeneric("ordering")
setGeneric("ordering",function(x,inv=FALSE)standardGeneric("ordering"))

setMethod("ordering","spam.chol.NgPeyton",function(x,inv=FALSE)
          {
            if (inv) return(x@invpivot) else return(x@pivot) })



setMethod("ordering","matrix",function(x,inv=FALSE)
          {
            if (dim(x)[1]!=dim(x)[2])
              stop("ordering is defined for square matrices only")
            if(inv)return(dim(x)[1]:1) else return(1:dim(x)[1]) })

setMethod("ordering","spam",function(x,inv=FALSE)
          {
            if (dim(x)[1]!=dim(x)[2])
              stop("ordering is defined for square matrices only")
            if(inv)return(dim(x)[1]:1) else return(1:dim(x)[1]) })



update.spam.chol.NgPeyton <- function(object,x,...){
  nrow <- object@dimension[1]
  if (!is.spam(x))
    stop("Covariance should be a 'spam' object.")
  if ((x@rowpointers[nrow+1]-1) != object@nnzA)
    stop("Updated covariance entries do not match length of original one.") 

  u <- .Fortran("updatefactor",
                nrow,
                object@nnzA,
                d =  as.double(x@entries),  jd = x@colindices,
                id = x@rowpointers,
                object@invpivot,   object@pivot,
                lindx=object@colindices,    xlindx=object@colpointers,
                nsuper=as.integer( length(object@supernodes)-1),
                entries = vector("double",length(object@entries)), #lnz
                rowpointers = object@rowpointers,#xlnz
                snode=object@snmember,  xsuper=object@supernodes,
                cachesize=object@memory[3],
                ierr = 0L,         
                NAOK = .Spam$NAOK,PACKAGE="spam")

  if(u$ierr>1) stop("Internal error in 'update.spam.chol.NgPeyton' code ", u$ierr,call.=FALSE)

  if(u$ierr == 1) {
    if (.Spam$cholupdatesingular == "null")
      return(NULL)
    else if (.Spam$cholupdatesingular == "error")
      stop("Singularity problem when updating a Cholesky Factor.")
    else if (.Spam$cholupdatesingular == "warning")
      warning("Singularity problem when updating a Cholesky Factor.\n'object' not updated.")
    else
      stop("'cholupdatesingular' should be 'error', 'null' or 'warning'.")
  }  else {
    slot(object, "entries", check = FALSE) <- u$entries
  }
  invisible(object)
}


chol.spam <- function(x, pivot = "MMD",
                      method="NgPeyton",
                      memory=list(),
                      eps = .Spam$eps, ...){

  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
  nrow <- x@dimension[1]
  nnzA <- as.integer( x@rowpointers[nrow+1]-1)
  if(nrow!=x@dimension[2]) stop("non-square matrix in 'chol'",call.=FALSE)

  if (any( diag.of.spam(x, nrow, nrow) < .Spam$eps))
    stop("Input matrix to 'chol' not positive definite (up to eps)",call.=FALSE)

  
  if(.Spam$cholsymmetrycheck) {
    test <- isSymmetric.spam(x, tol = eps*100) 
    if (!isTRUE(test))
      stop("Input matrix to 'chol' not symmetric (up to 100*eps)",call.=FALSE)
  }
  
  if (method != "NgPeyton")
    warning(gettextf("method = '%s' is not supported. Using 'NgPeyton'",
                     method), domain = NA)

  if (length(pivot)==1) {
    if (pivot==FALSE) {
      doperm <- 0L
      pivot <- seq_len(nrow)
    } else if(pivot==TRUE) {
      doperm <- 1L
      pivot <- vector("integer",nrow)
    } else {
      doperm <- as.integer( switch(match.arg(pivot,c("MMD","RCM")),MMD=1,RCM=2))
      pivot <- vector("integer",nrow)
    }
  } else  if (length(pivot)==nrow) {
    doperm <- 0L
    if (!is.integer(pivot[1]))
      pivot <- as.vector(pivot,"integer")
    if (.Spam$cholpivotcheck) {
      checkpivot(pivot,nrow)
    }
 } else stop("'pivot' should be 'MMD', 'RCM' or a valid permutation")



  ### IMPROVEME get better parameter values
  nnzcfact <- c(5,1,5)
  nnzRfact <- c(5,1,2)
  # nnzcolindices = length of array holding the colindices 
  if(is.null(memory$nnzcolindices))  {
    nnzcolindices <- ifelse((nnzA/nrow < 5), # very sparse matrix
                            max(1000,nnzA*(1.05*nnzA/nrow-3.8)),
                            nnzA)*nnzcfact[doperm+1]
    nnzcolindices <- max(nnzcolindices,nnzA)
 }else {
    nnzcolindices <- max(memory$nnzcolindices,nnzA)
    memory$nnzcolindices <- NULL
  }
  # nnzR = length of array holding the nonzero values of the factor 
  if(is.null(memory$nnzR))    nnzR <- min(max(4*nnzA,floor(.4*nnzA^1.2))*nnzRfact[doperm+1],nrow*(nrow+1)/2)  else {
    nnzR <- memory$nnzR
    memory$nnzR <- NULL
  }
  if(is.null(memory$cache))    cache <- 512  else {
    cache <- memory$cache 
    memory$cache <- NULL
  }

  if (length( memory)>0 )
    warning("The component(s) ", paste("'",names(memory),"'",sep='',collapse=","),
            " of the argument 'memory'\npassed to function 'chol' not meaningful and hence ignored.",call.=FALSE)
  nnzR <- as.integer(nnzR)
  nnzcolindices <- as.integer(nnzcolindices)
  z <- .Fortran("cholstepwise",
                nrow = nrow,nnzA = nnzA,
                d =  as.double(x@entries),jd = x@colindices,id = x@rowpointers,
                doperm = doperm,invp = vector("integer",nrow), perm = pivot,
                nnzlindx = vector("integer",1),             
                nnzcolindices = as.integer(nnzcolindices),
                lindx = vector("integer",nnzcolindices),     
                xlindx = vector("integer",nrow+1),     #
                nsuper = vector("integer",1),          #
                nnzR = as.integer(nnzR),#
                lnz = vector("double",nnzR),        #
                xlnz = vector("integer",nrow+1),     #
                snode = vector("integer",nrow),
                xsuper = vector("integer",nrow+1),   
                cachesize = as.integer(cache),
                ierr = 0L,          
                NAOK = .Spam$NAOK,PACKAGE="spam")

  if(z$ierr == 1) stop("Singularity problem when calculating the Cholesky factor.") 
  if(z$ierr == 6) stop("Inconsitency in the input",call.=FALSE)

  while( z$ierr>1) {
    if(z$ierr == 4) {
      tmp <- ceiling(nnzR*.Spam$cholincreasefactor[1])
      warning("Increased 'nnzR' with 'NgPeyton' method\n",
                    "(currently set to ",tmp," from ",nnzR,")",call.=FALSE)
      nnzR <- tmp
    }
    if(z$ierr == 5) {
      tmp <- ceiling(nnzcolindices*.Spam$cholincreasefactor[2])
      warning("Increased 'nnzcolindices' with 'NgPeyton' method\n",
         "(currently set to ",tmp," from ",nnzcolindices,")",call.=FALSE)
      nnzcolindices <- tmp
    }
    z <- .Fortran("cholstepwise",
                  nrow = nrow,nnzA = as.integer(x@rowpointers[nrow+1]-1),
                  d =  as.double(x@entries),jd = x@colindices,id = x@rowpointers,
                  doperm = doperm,invp = vector("integer",nrow), perm = pivot,
                  nnzlindx = vector("integer",1),             
                  nnzcolindices = as.integer(nnzcolindices),
                  lindx = vector("integer",nnzcolindices),     
                  xlindx = vector("integer",nrow+1),     #
                  nsuper = vector("integer",1),          #
                  nnzR = as.integer(nnzR),#
                  lnz = vector("double",nnzR),        #
                  xlnz = vector("integer",nrow+1),     #
                  snode = vector("integer",nrow),
                  xsuper = vector("integer",nrow+1),   
                  cachesize = as.integer(cache),
                  ierr = 0L,          
                  NAOK = .Spam$NAOK,PACKAGE="spam")
    
    if(z$ierr == 1) stop("Singularity problem when calculating the Cholesky factor.") 
  }
  nnzR <- as.integer(z$xlnz[length(z$xlnz)]-1)

  newx <- new("spam.chol.NgPeyton")
  slot(newx,"entries",check=FALSE) <- z$lnz[1:nnzR]
  slot(newx,"colindices",check=FALSE) <- z$lindx[1:z$nnzlindx]
  slot(newx,"colpointers",check=FALSE) <- z$xlindx[1:(z$nsuper+1)]   ########!!!!!!!
  slot(newx,"rowpointers",check=FALSE) <- z$xlnz
  slot(newx,"dimension",check=FALSE) <- c(nrow,nrow)
  slot(newx,"pivot",check=FALSE) <- z$perm
  slot(newx,"invpivot",check=FALSE) <- z$invp
  slot(newx,"supernodes",check=FALSE) <- z$xsuper[1:(z$nsuper+1)]
  slot(newx,"snmember",check=FALSE) <- z$snode
  slot(newx,"memory",check=FALSE) <- as.integer(c(nnzcolindices,z$nnzR,cache))
  slot(newx,"nnzA",check=FALSE) <- nnzA
  invisible(newx)
}

solve.spam <- function (a, b,  Rstruct = NULL, ...) {
  nrow <- a@dimension[1]
  ncol <- a@dimension[2]
  if (ncol != nrow)      stop("only square matrices can be inverted")

  if (missing(b)) {
    b <- diag(1, ncol)
  }  else {
    if(!is.matrix(b)) b <- as.matrix(b)
  }
  p <- dim(b)[2]
  if(nrow!=dim(b)[1])stop("'b' must be compatible with 'a'")
  
  # if we have a spam matrix, we calculate the Cholesky factor
  if (is(a,"spam"))
    if (is(Rstruct, "spam.chol.NgPeyton")) 
        a <- update.spam.chol.NgPeyton(Rstruct, a, ...)
    else a <- chol.spam(a, ...)


  if (is(a,"spam.chol.NgPeyton")) {
      # The following is a fast way to perform:
      #     z <- backsolve(a,forwardsolve( t(a),b))
    nsuper <- as.integer( length(a@supernodes)-1)
    z <- .Fortran("backsolves", m = nrow,
                  nsuper, p, a@colindices,
                  a@colpointers, as.double(a@entries),
                  a@rowpointers, a@invpivot, a@pivot,
                  a@supernodes, vector("double",nrow), sol = vector("double",nrow*p),
                  as.vector(b,"double"),
                  NAOK = .Spam$NAOK,PACKAGE = "spam")$sol
  } else z <- backsolve(a, forwardsolve( t(a),b))
    # see the helpfile for a comment about the 't(a)' construct.
  
  if ( p!=1)    dim(z) <- c(nrow,p)
  return( z)
}

chol2inv.spam <- function (x, ...) {
  nrow <- x@dimension[1]
  
  if (is(x,"spam.chol.NgPeyton")) {
    y <- vector("double",nrow*nrow)
    y[1L + 0L:(nrow - 1L) * (nrow + 1L)] <- 1.0

    z <- .Fortran("backsolves", m = nrow,
                  as.integer( length(x@supernodes)-1), nrow, x@colindices,
                  x@colpointers, as.double(x@entries),
                  x@rowpointers, x@invpivot, x@pivot,
                  x@supernodes, vector("double",nrow), sol = vector("double",nrow*nrow), y,
                  NAOK = .Spam$NAOK,PACKAGE = "spam")$sol
    dim(z) <- c(nrow,nrow)
  } else z <- backsolve.spam(x, forwardsolve.spam( t(x), diag(nrow)))
  return( z)
}

backsolve.spam <- function(r, x,...){#, k = NULL, upper.tri = NULL, transpose = NULL){
# r: spam.chol.NgPeyton structure as returned by chol.spam or a spam object
# x: rhs a vector or a matrix in dense form
# dimensions:  ( m x n) ( n x p) 
  m <- r@dimension[1]
  if(is.vector(x)) {
    n <- length(x)
    p <- 1L
  } else {
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
  }

  # we separate between "spam.chol.NgPeyton" and "spam"
  if (is(r,"spam.chol.NgPeyton")) {
    if (n!=m) stop("Cholesky factor 'r' not compatible with 'x'")
    nsuper <- as.integer( length(r@supernodes)-1)
    if (!.Spam$dopivoting) {
      z <- .Fortran("backsolve", m, nsuper, p, r@colindices,
                    r@colpointers, as.double(r@entries), r@rowpointers, 
                    r@supernodes, sol = vector("double",m*p),
                    NAOK = .Spam$NAOK,
                    PACKAGE="spam")$sol
    }else{
      z <- .Fortran("pivotbacksolve", m, nsuper, p, r@colindices,
                    r@colpointers,  as.double(r@entries),
                    r@rowpointers,   r@invpivot, r@pivot,
                    r@supernodes, vector("double",m),
                    sol = vector("double",m*p), as.double(x),
                    NAOK = .Spam$NAOK,PACKAGE="spam")$sol
    }
  } else {
    if (n!=m) stop("Triangular matrix 'r' not compatible with 'x'")
    # solve R sol = x
    z <- .Fortran("spamback",
                  m=m,p,sol = vector("double",m*p),x=as.vector(x,"double"),
                  al=as.double(r@entries),jal=r@colindices,
                  ial=r@rowpointers,
                  NAOK = .Spam$NAOK,PACKAGE="spam")
    if (z$m<0) stop(gettextf("singular matrix in 'backsolve'. Last zero in diagonal [%d]",
            -z$m), domain = NA)
     else z <- z$sol
   
  }
  
  if (p>1)     dim(z) <- c(m,p)
  return(z)
}

forwardsolve.spam <- function(l, x,...){#, k = NULL, upper.tri = NULL, transpose = NULL){
#  l: spam.chol.NgPeyton structure as returned by chol.spam
#         or an ordinary lower triangular spam matrix
#  x: rhs a vector a matrix in dense form
#  dimensions:  ( m x n) ( n x p)
#  if (!any(is.null(c(upper.tri,k,transpose ))))
#    warning("'k', 'upper.tri' and 'transpose' argument do not have any effect here")

  
  m <- l@dimension[1]
  if(is.vector(x)) {
    n <- length(x)
    p <- 1L
  } else {
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
  }

  # we separate between "spam.chol.NgPeyton" and "spam"
  if (is(l,"spam.chol.NgPeyton")) {
    if(n!=m) stop("Cholesky factor 'l' not compatible with 'x'")
    nsuper <- as.integer( length(l@supernodes)-1)
    if (!.Spam$dopivoting) {
      z <- .Fortran("forwardsolve", m, nsuper, p, l@colindices,
                    l@colpointers, as.double(l@entries), l@rowpointers, 
                    l@supernodes, sol = vector("double",m*p),
                    NAOK = .Spam$NAOK,PACKAGE="spam")$sol
    }else{
      z <- .Fortran("pivotforwardsolve", m, nsuper, p, l@colindices,
                    l@colpointers,  as.double(l@entries),
                    l@rowpointers,   l@invpivot, l@pivot,
                    l@supernodes, vector("double",m),
                    sol = vector("double",m*p), as.double(x),
                    NAOK = .Spam$NAOK,PACKAGE="spam")$sol
    }
  } else {
    if (n!=m) stop("Triangular matrix 'l' not compatible with 'x'")
    # solve L sol = x
    z <- .Fortran("spamforward",
                  m=m,p,sol = vector("double",m*p),x=as.vector(x,"double"),
                  al=as.double(l@entries),jal=l@colindices,
                  ial=l@rowpointers,
                  NAOK = .Spam$NAOK,PACKAGE="spam")
    if (z$m<0) stop(gettextf("singular matrix in 'forwardsolve'. First zero in diagonal [%d]",
            -z$m), domain = NA)
    else z <- z$sol
  }
  if (p>1)
    dim(z) <- c(m,p)
  return(z)
}


setMethod("chol","spam", chol.spam)
setMethod("solve","spam",solve.spam)
setMethod("chol2inv","spam", chol2inv.spam)
setMethod("chol2inv","spam.chol.NgPeyton", chol2inv.spam)

setMethod("backsolve","spam",#signature(r="spam",x='ANY'),
          backsolve.spam)
setMethod("backsolve","spam.chol.NgPeyton",#signature(r="spam.chol.NgPeyton",x='ANY'),
          backsolve.spam,sealed=TRUE)
#setMethod("backsolve","spam.chol.NgPeyton",    backsolve.spam)
setMethod("forwardsolve","spam",               forwardsolve.spam)
setMethod("forwardsolve","spam.chol.NgPeyton", forwardsolve.spam)

######################################################################
######################################################################

determinant.spam <- function(x, logarithm = TRUE, pivot = "MMD",method="NgPeyton",
                              memory=list(),eps = .Spam$eps, ...){
  
  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
  logdet <- list()
 #### start from above 
  nrow <- x@dimension[1]
  nnzA <- as.integer( x@rowpointers[nrow+1]-1)
  if(nrow!=x@dimension[2]) stop("non-square matrix in 'chol'",call.=FALSE)

  if(.Spam$cholsymmetrycheck) {
    test <- isSymmetric.spam(x, tol = eps*100) 
    if (!isTRUE(test))
      stop("Input matrix to 'chol' not symmetric (up to 100*eps)",call.=FALSE)
  }
  
  if (method != "NgPeyton")
    warning(gettextf("method = '%s' is not supported. Using 'NgPeyton'",
                     method), domain = NA)

  if (length(pivot)==nrow) {
    doperm <- 0L
    pivot <- as.vector(pivot,"integer")
    if (.Spam$cholpivotcheck) {
      checkpivot(pivot,nrow)
    }
  } else if (length(pivot)==1) {
    if (pivot==FALSE) {
      doperm <- 0L
      pivot <- seq_len(nrow)
    } else if(pivot==TRUE) {
      doperm <- 1L
      pivot <- vector("integer",nrow)
    } else {
      doperm <- as.integer( switch(match.arg(pivot,c("MMD","RCM")),MMD=1,RCM=2))
      pivot <- vector("integer",nrow)
    }
  } else stop("'pivot' should be 'MMD', 'RCM' or a permutation")


  ### IMPROVEME get better parameter values
  nnzcfact <- c(5,1,5)
  nnzRfact <- c(5,1,2)
  # nnzcolindices = length of array holding the colindices 
  if(is.null(memory$nnzcolindices))  {
    nnzcolindices <- ifelse((nnzA/nrow < 5), # very sparse matrix
                            max(1000,nnzA*(1.05*nnzA/nrow-3.8)),
                            nnzA)*nnzcfact[doperm+1]
    nnzcolindices <- max(nnzcolindices,nnzA)
 }else {
    nnzcolindices <- max(memory$nnzcolindices,nnzA)
    memory$nnzcolindices <- NULL
  }
  # nnzR = length of array holding the nonzero values of the factor 
  if(is.null(memory$nnzR))    nnzR <- min(max(4*nnzA,floor(.4*nnzA^1.2))*nnzRfact[doperm+1],nrow*(nrow+1)/2)  else {
    nnzR <- memory$nnzR
    memory$nnzR <- NULL
  }
  if(is.null(memory$cache))    cache <- 64  else {
    cache <- memory$cache 
    memory$cache <- NULL
  }

  if (length( memory)>0 )
    warning("The component(s) ", paste("'",names(memory),"'",sep='',collapse=","),
            " of the argument 'memory'\npassed to function 'chol' not meaningful and hence ignored.",call.=FALSE)
  
  z <- .Fortran("cholstepwise",
                nrow = nrow,nnzA = as.integer(x@rowpointers[nrow+1]-1),
                d =  as.double(x@entries),jd = x@colindices,id = x@rowpointers,
                doperm = doperm,invp = vector("integer",nrow), perm = pivot,
                nnzlindx = vector("integer",1),             
                nnzcolindices = as.integer(nnzcolindices),
                lindx = vector("integer",nnzcolindices),     
                xlindx = vector("integer",nrow+1),     #
                nsuper = vector("integer",1),          #
                nnzR = as.integer(nnzR),#
                lnz = vector("double",nnzR),        #
                xlnz = vector("integer",nrow+1),     #
                snode = vector("integer",nrow),
                xsuper = vector("integer",nrow+1),   
                cachesize = as.integer(cache),
                ierr = 0L,          
                NAOK = .Spam$NAOK, PACKAGE = "spam")

  if(z$ierr == 1) stop("Singularity problem when calculating the Cholesky factor.") 
  if(z$ierr == 6) stop("Inconsitency in the input",call.=FALSE)

  while( z$ierr>1) {
    if(z$ierr == 4) {
      warning("Increased 'nnzR' with 'NgPeyton' method\n",
              "(currently set to ",nnzR," from ",ceiling(nnzR*.Spam$cholpar[1]),")",call.=FALSE)
      nnzR <- ceiling(nnzR*.Spam$nnzRinc)
    }
    if(z$ierr == 5) {
      warning("Increased 'nnzcolindices' with 'NgPeyton' method\n",
         "(currently set to ",nnzcolindices," from ",ceiling(nnzcolindices*.Spam$cholpar[2]),")",call.=FALSE)
      nnzcolindices <- ceiling(nnzcolindices*.Spam$cholpar[2])
    }
    z <- .Fortran("cholstepwise",
                  nrow = nrow,nnzA = as.integer(x@rowpointers[nrow+1]-1),
                  d =  as.double(x@entries),jd = x@colindices,id = x@rowpointers,
                  doperm = doperm,invp = vector("integer",nrow), perm = pivot,
                  nnzlindx = vector("integer",1),             
                  nnzcolindices = as.integer(nnzcolindices),
                  lindx = vector("integer",nnzcolindices),     
                  xlindx = vector("integer",nrow+1),     #
                  nsuper = vector("integer",1),          #
                  nnzR = as.integer(nnzR),#
                  lnz = vector("double",nnzR),        #
                  xlnz = vector("integer",nrow+1),     #
                  snode = vector("integer",nrow),
                  xsuper = vector("integer",nrow+1),   
                  cachesize = as.integer(cache),
                  ierr = 0L,          
                  NAOK = .Spam$NAOK, PACKAGE = "spam")
    
  }
 #### end from above 
  if(z$ierr == 1) {
                                        # all other errors trapped 
      warning("singularity problem or matrix not positive definite",call.=FALSE)
      logdet$modulus <- NA
   } else{
    tmp <- 2* sum( log( z$lnz[ z$xlnz[ -(z$nrow+1)]]))
    if (logarithm) logdet$modulus <- tmp else logdet$modulus <- exp(tmp)
  }

  attr(logdet$modulus,"logarithm") <- logarithm
  
  logdet$sign <- ifelse(z$ierr == 1,NA,1)
  attr(logdet,"class") <- "det"
  
  return(logdet)
}

determinant.spam.chol.NgPeyton <- function(x, logarithm = TRUE,...)
{
  logdet <- list()

  
  tmp <- sum( log(x@entries[ x@rowpointers[-(x@dimension[1]+1)]]))
  if (logarithm) logdet$modulus <- tmp else logdet$modulus <- exp(tmp)
 
  attr(logdet$modulus,"logarithm") <- logarithm
  
  logdet$sign <- 1
  attr(logdet,"class") <- "det"
  
  return(logdet)
}


setMethod("determinant","spam",               determinant.spam)
setMethod("determinant","spam.chol.NgPeyton", determinant.spam.chol.NgPeyton)

######################################################################
########################################################################

    
"as.matrix.spam.chol.NgPeyton" <- function(x,...){
  nrow <- x@dimension[1]
  nnzR <- x@rowpointers[nrow+1]-1
  newx <- new("spam")
  nsuper <- as.integer( length(x@supernodes)-1)
  xcolindices <- .Fortran('calcja',
                           nrow, nsuper, x@supernodes, x@colindices, x@colpointers, x@rowpointers,
                           xja=vector("integer",nnzR),
                           NAOK = .Spam$NAOK,PACKAGE = "spam")$xja
  return(array(.Fortran("spamcsrdns",
                 nrow=nrow,
                 entries=as.double(x@entries),
                 colindices=xcolindices,
                 rowpointers=x@rowpointers,
                 res=vector("double",nrow*nrow),  
                 NAOK=.Spam$NAOK,PACKAGE = "spam")$res,
               c(nrow,nrow))      # we preserve dimensions
         )
}



setMethod("as.matrix","spam.chol.NgPeyton",as.matrix.spam.chol.NgPeyton)
setMethod("as.vector","spam.chol.NgPeyton",
          function(x){
            as.vector.spam(as.spam.chol.NgPeyton(x))
          })


########################################################################
#  force to spam matrices. Would not be required with inheritance

setMethod("image","spam.chol.NgPeyton",
          function(x,cex=NULL,...){
            image.spam(as.spam.chol.NgPeyton(x),cex=cex,...)
          })


setMethod("display","spam.chol.NgPeyton",
          function(x,...){
            display.spam(as.spam.chol.NgPeyton(x),...)
          })

setMethod("t","spam.chol.NgPeyton",
          function(x){
            t.spam(as.spam.chol.NgPeyton(x))
          })

setMethod("chol","spam.chol.NgPeyton",
          function(x){
           x
          })
########################################################################

### system.time({ for (i in 1:1000) x=1:1000000})           # 8.820 
### system.time({ for (i in 1:1000) x=seq(length=1000000)}) # 8.397
### system.time({ for (i in 1:1000) x=seq_len(1000000)})    # 8.628
### system.time({ for (i in 1:1000) x=seq.int(1000000)})    # 8.944

### system.time({ for (i in 1:100000) x=1:10000})           # 2.161
### system.time({ for (i in 1:100000) x=seq(length=10000)}) # 3.288
### system.time({ for (i in 1:100000) x=seq_len(10000)})    # 2.060
### system.time({ for (i in 1:100000) x=seq.int(10000)})    # 2.249
