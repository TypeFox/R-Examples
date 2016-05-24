setClass("mhp",  # "mhp" ==  "Multivariate HyperParameters"
         representation = representation(
           M      = "matrix",
           B      = "array",
           levels = "character",
           names  = "character"
           )
         )

"mhp" <- function(M,B,levels=NULL,names=NULL){
  if(is.null(names)){
      names <- dimnames(B)[[1]]
  }
  
  if(is.null(levels)){
    if(is.null(rownames(M))){
      levels <- dimnames(B)[[3]]
    } else {
      levels <- rownames(M)
    }
  }

  dimnames(M) <- NULL
  dimnames(B) <- NULL  # these 2 lines needed to ensure the call to
                       # new() passes .mhp_valid()
  
  new("mhp", M=M, B=B, levels=levels, names=names)  # This is the
                                                    # *ONLY*
                                                    # occurrence of
                                                    # new("mhp",...)
}

"is.mhp" <- function(x){is(x,"mhp")}

setGeneric("levels",function(x){standardGeneric("levels")})
setMethod("levels","mhp",function(x){x@levels})

"M" <- function(x){
  stopifnot(is.mhp(x))
  out <-   x@M
  rownames(out) <- levels(x)
  colnames(out) <- levels(x)
  return(out)
} 
          
"B" <- function(x){
  stopifnot(is.mhp(x))
  out <- x@B
  dimnames(out) <- list(names(x),names(x),levels(x))
  return(out)
} 


setGeneric("names")
setMethod("names","mhp",function(x){x@names})

".mhp_valid" <- function(object){
  M <- object@M
  B <- object@B  # NB: not 'B <- B(object)' because this puts in
                 # dimnames, which should be NULL (to avoid possible
                 # conflict between: dimnames(B)[[3]] and the 'types'
                 # argument; and rownames(M) and the 'names' argument.

  levs <- levels(object)
  
  if(!ipd(M)){
    return("M not positive definite")
  } else if(!is.null(dimnames(M))){
    return("M has names (supply these through the 'names' argument)")
  } else if(nrow(M) != length(levs)){
    return("size of M incompatible with length of names")
  } else if(!all(apply(B,3,ipd))){
    return("non-positive-definite slice of B")
  } else if(!is.null(dimnames(B))){
    return("B has dimnames (supply these through the 'names' argument)")
  } else if(dim(B)[3] != length(levs)){
    return("size of B incompatible with length of 'levels' argument")
  } else if(dim(B)[1] != length(object@names)){
    return("size of B incompatible with length of 'names' argument")
  } else {
    return(TRUE)
  }
}

setValidity("mhp", .mhp_valid)

setGeneric("types",function(x){standardGeneric("types")})
setMethod("types","mhp",function(x){stop("do not use types() on an mhp object.  You probably mean levels(). Function types() only works for mdm objects.")})

setGeneric("types<-",function(x,value){standardGeneric("types<-")})
setGeneric("names<-")

"M<-" <- function(x,value){
  stopifnot(is.mhp(x))
  mhp(M=value, B=B(x), levels(x), names=names(x))
} 

"B<-" <- function(x,value){
  stopifnot(is.mhp(x))
  mhp(M=M(x), B=value, levels(x), names=names(x))
}

setMethod("types<-","mhp",function(x,value){
  stop("do not use types<-() on an mhp object.  You probably mean levels<-(). Function types<-() only works for mdm objects.")
} )

setMethod("names<-","mhp",function(x,value){
  mhp(M=M(x), B=B(x), levels=levels(x), names=value)
} )

".mhp_print" <- function(x){
  list(M=M(x),B=B(x))
}
    
"print.mhp" <- function(x, ...){
  jj <- .mhp_print(x, ...)
  print(jj)
  return(invisible(jj))
}

setMethod("show", "mhp", function(object){print.mhp(object)})

#following lines adapted from the Matrix package
setMethod("summary", signature(object = "mhp"),
          function(object, ...) {
            jj <- B(object)
            # set on-diagonal elements to zero:
            jj[as.matrix(expand.grid(lapply(as.list(dim(jj)[-1]),seq_len))[,c(1,1,2)])] <- 0
            if(all(jj==0)){
              diag <- TRUE   #ie all offdiagonal elements are zero
            } else {
              diag <- FALSE
            }
            if(!diag){
              return(object)
            } else {
              jj <- B(object)
              out <- list(
                          M = M(object),
                          sapply(seq_len(dim(jj)[3]), function(i){diag(jj[,,i])},simplify=FALSE)
                          )
              class(out) <- 'mhpSummary'
              return(out)
            }
          })

print.mhpSummary <- function (x, ...) {
  cat('overall covariance matrix M:\n\n')
  print(x[[1]])
  cat('\n\n')
  cat('each B matrix is diagonal, with entries\n\n')
  jj <- x[[-1]]
  names(jj) <- rownames(x[[1]])
  print(jj)
  return(invisible(x))
}
