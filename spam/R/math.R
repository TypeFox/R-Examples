# This is file ../spam/R/math.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     

# `"Ops"':
#      `"+"', `"-"', `"*"', `"/"', `"^"', `"%%"', `"%/%"'
#      `"&"', `"|"', `"!"'
#      `"=="', `"!="', `"<"', `"<="', `">="', `">"'


#     `Math' `"abs"', `"sign"', `"sqrt"', `"ceiling"', `"floor"',
#          `"trunc"', `"cummax"', `"cummin"', `"cumprod"', `"cumsum"',
#          `"log"', `"log10"', `"log2"', `"log1p"', `"acos"', `"acosh"',
#          `"asin"', `"asinh"', `"atan"', `"atanh"', `"exp"', `"expm1"',
#          `"cos"', `"cosh"', `"cospi"', `"sin"', `"sinh"', `"sinpi"',
#          `"tan"', `"tanh"', `"tanpi"', `"gamma"', `"lgamma"',
#          `"digamma"', `"trigamma"'

#     `Math2' `"round"', `"signif"'

#     `Summary' `"max"', `"min"', `"range"', `"prod"', `"sum"', `"any"', `"all"'


##############
# Unary operators "+", "-" and "!" are handled with e2 missing...
#
# Currently, "+", "-" are handled...
setMethod("!",signature(x="spam"),    function(x){
    if(.Spam$structurebased) {
        x@entries <- as.double(callGeneric(x@entries))
        x
    } else {
        inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)), prod(dim(x)))
        spam(as.double( callGeneric(as.matrix(x))), nrow=nrow(x))
    }
})
          
setMethod("+",signature(e1="spam",e2="missing"), function(e1) e1 )                
setMethod("-",signature(e1="spam",e2="missing"), function(e1) { e1@entries <- -e1@entries; e1} )                
              
#     `Math2' :

setMethod("Math2",signature(x = "spam", digits = "ANY"),
          function(x, digits){ x@entries <- callGeneric(x@entries, digits = digits); x })

#     `Math' :

setMethod("Math","spam", function(x){
    if(.Spam$structurebased) {
        x@entries <- callGeneric(x@entries)
        x
    }else{
        x@entries <- callGeneric(x@entries)
        as.spam.spam( x)  
    }
})

#     `Math', where we pass to matrix first...

spam_Math <- function(x) {
    if(.Spam$structurebased) {
        x@entries <- callGeneric(x@entries)
        x
    }else{
        inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)), prod(dim(x)))
        as.spam(callGeneric(as.matrix(x)))
    }}


setMethod("exp","spam", spam_Math )
setMethod("log10","spam", spam_Math )
setMethod("log2","spam", spam_Math )
# from ?log: Do not set S4 methods on `logb' itself.
# special case to set base...          
setMethod("log","spam", function(x,...) {
    if(.Spam$structurebased) {
        x@entries <- callGeneric(x@entries,...)
        x
    }else{
        inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)), prod(dim(x)))
        as.spam(callGeneric(as.matrix(x),...))
    }}
          )


          
setMethod("cos","spam", spam_Math )
#setMethod("cospi","spam", spam_Math )
setMethod("cosh","spam", spam_Math )
setMethod("acosh","spam", spam_Math )
setMethod("acos","spam", spam_Math )

setMethod("gamma","spam", spam_Math )
setMethod("digamma","spam", spam_Math )
setMethod("trigamma","spam", spam_Math )
setMethod("lgamma","spam", spam_Math )

setMethod("cummax","spam", spam_Math )
setMethod("cummin","spam", spam_Math )
setMethod("cumprod","spam", spam_Math )
setMethod("cumsum","spam", spam_Math )


#     `Summary' :
setMethod("Summary","spam", function(x,...,na.rm=FALSE){
    if(.Spam$structurebased) {
        callGeneric(x@entries,...,na.rm=na.rm) 
    }else{
        if ( prod( x@dimension) == length( x@entries)) {
            callGeneric(x@entries,...,na.rm=na.rm) 
        } else {
            callGeneric(c(0,x@entries),...,na.rm=na.rm)
        }
    }
}
          )          

logical_Summary <- function( x,...,na.rm=FALSE){
    if(.Spam$structurebased) {
        callGeneric(as.logical(x@entries),...,na.rm=na.rm) 
    }else{
        if ( prod( x@dimension) == length( x@entries)) {
            callGeneric(as.logical(x@entries),...,na.rm=na.rm) 
        } else {
            callGeneric(as.logical(c(0,x@entries)),...,na.rm=na.rm)
        }
    }
}
      

setMethod("any","spam", logical_Summary)
setMethod("all","spam", logical_Summary)


################################################################################################################################################################################################################################################################################################
#     `Ops' `"Arith"', `"Compare"', `"Logic"'



#     `Logic' `"&"', `"|"'.
        
                                       
"spam_Logic_vectorspam" <- function(e1, e2) {
    if(.Spam$structurebased) {
        if(identical(length(e1),1L) | identical(length(e1), length(e2@entries))) {
            e2@entries <- as.double( callGeneric(e1, e2@entries))
            return(e2)
        }
        if( length(e1) == prod(e2@dimension))
            return( as.spam( callGeneric(e1, as.matrix(e2))) )
        stop(gettextf("incompatible lengths for %s operation.", sQuote(.Generic)))
    }  else {
        inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)), prod(dim(e2)))
        return( as.spam( callGeneric(e1, as.matrix(e2))) )
    }
}

"spam_Logic_spamvector" <- function(e1, e2)  {
    if(.Spam$structurebased) {
        if(identical(length(e2),1L) | identical(length(e2), length(e1@entries))) {
            e1@entries <- as.double( callGeneric(e1@entries, e2))
            return(e1)
        }
        if( length(e2)== prod(e1@dimension))
            return( as.spam( callGeneric(as.matrix(e1), e2)) )
        stop(gettextf("incompatible lengths for %s operation.", sQuote(.Generic)))
    }  else {
        inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)), prod(dim(e1)))
        return( as.spam( callGeneric(as.matrix(e1), e2)) )
    }
}
        

setMethod("|",signature(e1="spam",e2="spam"), 
          function(e1,e2){ z <- spam_add(e1,e2);z@entries <- rep(1,length(z@colindices));z})

setMethod("&",signature(e1="spam",e2="spam"), 
          function(e1,e2){ z <- spam_mult(e1,e2); z@entries <- rep(1,length(z@colindices));z})
setMethod("Logic",signature(e1="spam",e2="vector"), spam_Logic_spamvector)
setMethod("Logic",signature(e1="vector",e2="spam"), spam_Logic_vectorspam)

##################################################################################################
#     `Compare' `"=="', `">"', `"<"', `"!="', `"<="', `">="'                                     
"spam_Compare" <- function(e1,e2) {
    inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)),  max(prod(dim(e1)), prod(dim(e2))))
    as.spam( callGeneric( as.matrix(e1), as.matrix(e2))  )            
}

"spam_Compare_spamvector" <- function(e1, e2){
    if(.Spam$structurebased) {
        if(identical(length(e2),1L) | identical(length(e2), length(e1@entries))) {
            e1@entries <- as.double(callGeneric(e1@entries, e2))
            return(e1)
        }
        if( length(e2)== prod(e1@dimension))
            return( as.spam( callGeneric(as.matrix(e1), e2)) )
        stop(gettextf("incompatible lengths for %s operation.", sQuote(.Generic)))
    }  else {
        inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)),   prod(dim(e1)))
        return( as.spam( callGeneric(as.matrix(e1), e2)) )
    }
}
"spam_Compare_vectorspam" <- function(e1, e2) {
    if(.Spam$structurebased) {
        if(identical(length(e1),1L) | identical(length(e1), length(e2@entries))) {
            e2@entries <- as.double( callGeneric(e1, e2@entries))
            return(e2)
        }
        if( length(e1) == prod(e2@dimension))
            return( as.spam( callGeneric(e1, as.matrix(e2))) )
        stop(gettextf("incompatible lengths for %s operation.", sQuote(.Generic)))
     }  else {
        inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)),   prod(dim(e2)))
        return( as.spam( callGeneric(e1, as.matrix(e2))) )
    }
}


setMethod("Compare",signature(e1="spam",e2="spam"),   spam_Compare )
setMethod("Compare",signature(e1="spam",e2="vector"), spam_Compare_spamvector )
setMethod("Compare",signature(e1="vector",e2="spam"), spam_Compare_vectorspam )
##################################################################################################
#     `Arith': `"+"', `"-"', `"*"', `"^"', `"%%"', `"%/%"', `"/"'

"spam_Arith_vectorspam" <- function(e1, e2){    
#    cat("spam_Arith_vectorspam")
    if(.Spam$structurebased) {
        if(identical(length(e1),1L) | identical(length(e1), length(e2@entries))) {
            e2@entries <- callGeneric(e1, e2@entries)
            return(e2)
        }
        if( length(e1) == prod(e2@dimension))
            return( as.spam( callGeneric(e1, as.matrix(e2))) )
        stop(gettextf("incompatible lengths for %s operation.", sQuote(.Generic)))
    }  else {
        inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)), prod(dim(e1)))
         return( as.spam( callGeneric(e1, as.matrix(e2))) )
    }
}
"spam_Arith_spamvector" <- function(e1, e2){
#    cat("spam_Arith_spamvector")
    if(.Spam$structurebased) {
        if(identical(length(e2),1L) | identical(length(e2), length(e1@entries))) {
            e1@entries <- callGeneric(e1@entries, e2)
            return(e1)
        }
        if( length(e2)== prod(e1@dimension))
            return( as.spam( callGeneric(as.matrix(e1), e2)) )
        stop(gettextf("incompatible lengths for %s operation.", sQuote(.Generic)))
     }  else {
        inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)), prod(dim(e1)))
        return( as.spam( callGeneric(as.matrix(e1), e2)) )
    }
}
spam_Arith <- function(e1,e2) {
    inefficiencywarning( gettextf("This %s operation may be inefficient",sQuote(.Generic)),  max(prod(dim(e1)), prod(dim(e2))))
    as.spam( callGeneric( as.matrix(e1), as.matrix(e2)))
}
        
 

setMethod("Arith",signature(e1="spam",e2="spam"),   spam_Arith )
setMethod("Arith",signature(e1="spam",e2="vector"), spam_Arith_spamvector)
setMethod("Arith",signature(e1="vector",e2="spam"), spam_Arith_vectorspam)

setMethod("/",signature(e1="spam",e2="spam"), function(e1,e2){ "/"(e1,as.matrix(e2)) } )
setMethod("^",signature(e1="spam",e2="spam"), function(e1,e2){ "^"(e1,as.matrix(e2)) } )


######################################################################
# nz <- 128; ln <- nz^2; A <- spam(0,ln,ln); is <- sample(ln,nz); js <- sample(ln,nz);A[cbind(is,js)] <- 1:nz
# nz <- 128; ln <- nz^2; A <- spam(0,ln,ln); is <- sample(ln,ln); js <- sample(ln,ln);A[cbind(is,js)] <- 1:ln
# system.time(   spam:::.spam.addsparsefull(A,1)) ; system.time(   as.matrix.spam(A)+1.5)


#######################################################################
"spam_add" <- function(A, B, s=1)
{

    # cat("spam_add")
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
  if(ncol != B@dimension[2] || nrow != B@dimension[1])
    stop("non-conformable matrices")
   
  nzmax <- .Fortran("aplbdg",
                nrow,                ncol,
                A@colindices,               A@rowpointers,
                B@colindices,               B@rowpointers,
                vector("integer",nrow),nnz=vector("integer",1),vector("integer",ncol),
                NAOK=.Spam$NAOK,PACKAGE = "spam")$nnz

  z <- .Fortran("aplsb1",
                nrow,
                ncol,
                as.double(A@entries), A@colindices,  A@rowpointers,
                as.double(s),
                as.double(B@entries), B@colindices,  B@rowpointers,
                entries     = vector("double",nzmax),
                colindices  = vector("integer",nzmax),
                rowpointers = vector("integer",nrow+1),
                as.integer(nzmax+1),  ierr = vector("integer",1),
                NAOK=.Spam$NAOK,PACKAGE = "spam")
  if(z$ierr != 0) stop("insufficient space for sparse matrix addition")
  nz <- z$rowpointers[nrow+1]-1
  newz <- new("spam")
  slot(newz,"entries", check=FALSE) <- z$entries[1:nz]
  slot(newz,"colindices", check=FALSE) <- z$colindices[1:nz]
  slot(newz,"rowpointers", check=FALSE) <- z$rowpointers
  slot(newz,"dimension", check=FALSE) <- c(nrow,ncol)
  return(newz)
}

setMethod("+",signature(e1="spam",e2="spam"),  function(e1,e2){ spam_add(e1, e2)    })
setMethod("-",signature(e1="spam",e2="spam"),  function(e1,e2){ spam_add(e1, e2, -1)})



###############################################################################

"spam_mult" <- function(e1,e2)
{
#  if(is.vector(e1)) {
#    if(length(e1) == 1){
#      if(e1==0) return( spam(0,nrow(e2),ncol(e2)))
#      else{  # just a scalar
#        e2@entries <- e1*e2@entries
#        return(e2)
#      }
#    }  else if(length(e1) == nrow(e2))
#      return(diag.spam(e1) %*% e2)
#    else # length(e1) == ncol(e2) is not required
#      stop("e1 and e2 not conformable for efficient element-by-element multiplication")
#  }
#  else if(is.vector(e2)) {
#    if(length(e2) == 1){
#      if(e2==0)   return( spam(0,nrow(e1),ncol(e1)))
#      else {
#        e1@entries <- e2*e1@entries
#        return(e1)
#      }
#    }
#    else if(length(e2) == nrow(e1))
#      return(diag.spam(e2) %*% e1)
#    else
#      stop("e1 and e2 not conformable for efficient element-by-element multiplication")
#  }
#  if(is.matrix(e1))
#    e1 <- as.spam(e1)
#  else if(is.matrix(e2))
#    e2 <- as.spam(e2)
#  if(!(is.spam(e1) && is.spam(e2)))
#    stop("Arguments must be of class:  vector, matrix or spam")
  
  e1row <- e1@dimension[1]
  e1col <- e1@dimension[2]
  if(e1col != e2@dimension[2] | e1row != e2@dimension[1])
    stop("non-conformable matrices")
  nnzmax <- length(intersect(e1@colindices+e1col*(rep(1:e1row,diff(e1@rowpointers))-1),
                             e2@colindices+e2@dimension[2]*(rep(1:e2@dimension[1],diff(e2@rowpointers))-1)))+1
  z <- .Fortran("aemub",
                e1row,
                e1col,
                as.double(e1@entries), e1@colindices,  e1@rowpointers,
                as.double(e2@entries), e2@colindices,  e2@rowpointers,
                entries     = vector("double",nnzmax),
                colindices  = vector("integer",nnzmax),
                rowpointers = vector("integer",e1row+1),
                integer(e1col),
                double(e1col),
                as.integer(nnzmax),  ierr = vector("integer",1),
                NAOK=.Spam$NAOK,  PACKAGE = "spam")
  if(z$ierr != 0)      stop("insufficient space for element-wise sparse matrix multiplication")
  nnz <- z$rowpointers[e1row+1]-1
  if(identical(z$entries,0L)){#trap zero matrix
    z$colindices <- 1L
    z$rowpointers <- c(1L,rep(2L,e1row))
  }
  return(new("spam",entries=z$entries[1:nnz],colindices=z$colindices[1:nnz],rowpointers=z$rowpointers,
             dimension=c(e1row,e1col)))
}


setMethod("*",signature(e1="spam",e2="spam"), spam_mult)


##########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
"matrix_add_spammatrix" <- function(A,B){
#    cat("matrix_add_spammatrix")
  # A is sparse, B is full
                 #  if (missing(B)) return(A)
                 #  if (!is.numeric(B)) stop("numeric argument expected")
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
  pdim <- prod(nrow,ncol)
#  if (is.matrix(B)) {
    if(ncol != dim(B)[2] || nrow != dim(B)[1])
      stop("non-conformable matrices")
#  } else {
#    if(pdim%%length(B)!=0) {
#      stop("longer object length
#        is not a multiple of shorter object length")
#    } else  B <- rep(B,pdim %/% length(B))
#  }
  return(array( .Fortran("addsparsefull",
                          nrow,as.double(A@entries),A@colindices,
                          A@rowpointers,b=as.double(B),NAOK=.Spam$NAOK,PACKAGE = "spam"
                          )$b,c(nrow,ncol)))
 }

"matrix_sub_spammatrix" <- function(A,B){
#    cat("matrix_sub_spammatrix")
  # A is sparse, B is full
#  if (missing(B)) {
#    A@entries <- -A@entries
#    return(A)
#  }
#  if (!is.numeric(B)) stop("numeric argument expected")
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
#  pdim <- prod(nrow,ncol)
#  if (is.matrix(B)) {
    if(ncol != dim(B)[2] || nrow != dim(B)[1])
      stop("non-conformable matrices")
#  } else {
#    if(pdim %% length(B)!=0) {
#      stop("longer object length
#        is not a multiple of shorter object length")
#    } else  B <- rep(B,pdim %/% length(B))
#  }
#  if (!is.double(B[1]))  B <- as.double(B)
  return(array( .Fortran("subfullsparse",
                          nrow,ncol,as.double(A@entries),A@colindices,
                          A@rowpointers,b=as.double(B),NAOK=.Spam$NAOK,PACKAGE = "spam"
                          )$b,c(nrow,ncol)))
}

"matrix_sub_matrixspam" <- function(B,A){
#    cat("matrix_sub_spammatrix")
  # A is sparse, B is full
  if (!is.numeric(B)) stop("numeric argument expected")
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
#  pdim <- prod(nrow,ncol)
#  if (is.matrix(B)) {
    if(ncol != dim(B)[2] || nrow != dim(B)[1])
      stop("non-conformable matrices")
#  } else {
#    if(pdim %% length(B)!=0) {
#      stop("longer object length
#        is not a multiple of shorter object length")
#    } else  B <- rep(B,pdim %/% length(B))
#  }
#  if (!is.double(B[1]))  B <- as.double(B)
  return(array( .Fortran("subsparsefull",
                          nrow,as.double(A@entries),A@colindices,
                          A@rowpointers,b=as.double(B),NAOK=.Spam$NAOK,PACKAGE = "spam"
                          )$b,c(nrow,ncol)))
}
#
setMethod("+",signature(e1="spam",   e2="matrix"), function(e1,e2){ matrix_add_spammatrix(e1,e2)})
setMethod("+",signature(e1="matrix",   e2="spam"), function(e1,e2){ matrix_add_spammatrix(e2,e1)})
setMethod("-",signature(e1="matrix",   e2="spam"), function(e1,e2){ matrix_sub_matrixspam(e1,e2)})
setMethod("-",signature(e1="spam",   e2="matrix"), function(e1,e2){ matrix_sub_spammatrix(e1,e2)})

                                      
#"spam_division" <- function(e1,e2) { # Element-wise matrix division of two spams
#  if(is.numeric(e1) && length(e1) == 1)
#  { e2@entries <- e1/e2@entries
#    return(e2)
# } else if(is.numeric(e2) && length(e2) == 1) {
#    e1@entries <- e1@entries/e2
#    return(e1)
#  }
#  else if(is.spam(e1) || is.spam(e2) || is.matrix(e1) || is.matrix(e2)){
#        if(is.matrix(e1)) e1 <- as.spam(e1)
#        if(is.matrix(e2)) e2 <- as.spam(e2)
#        nrow <- e1@dimension[1]
#        ncol <- e1@dimension[2]
#        if(ncol != e2@dimension[2] | nrow != e2@dimension[1])
#          stop("matrices not conformable for element-by-element division")
#	nzmax <- length(unique(c(e1@colindices+ncol*(rep(1:nrow,diff(e1@rowpointers))-1),
#                                 e2@colindices+e2@dimension[2]*(rep(1:e2@dimension[1],diff(e2@rowpointers))-1))))+1
#        z <- .Fortran("_aedib_",   # does not order the colindicies upon return!
#                      nrow,
#                      ncol,
#                      as.integer(1),
#                      as.double(e1@entries), e1@colindices,  e1@rowpointers,
#                      as.double(e2@entries), e2@colindices,  e2@rowpointers,
#                      entries     = vector("double",nzmax),
#                      colindices  = vector("integer",nzmax),
#                      rowpointers = vector("integer",nrow+1),
#                      as.integer(nzmax),
#                      integer(ncol),
#                      double(ncol),
#                      ierr = vector("integer",1),
#                      NAOK=.Spam$NAOK,PACKAGE = "spam"
#                      )
#        if(z$ierr != 0) stop("insufficient space for element-wise sparse matrix division")
#        nz <- z$rowpointers[nrow+1]-1
#        return(new("spam",entries=z$entries[1:nz],colindices=z$colindices[1:nz],rowpointers=z$rowpointers,
#                   dimension=c(nrow,ncol)))
#    }
 

   
##"spam_exponent" <- function(e1, e2)
#{
#    nrow <- e1@dimension[1]
#    ncol <- e1@dimension[2]
#    if(ncol != e2@dimension[2] | nrow != e2@dimension[1])
#        stop("matrices not conformable for element-wise exponentiation ")
#    nzmax <- length(unique(c(e1@colindices+ncol*(rep(1:nrow,diff(e1@rowpointers))-1),
#                             e2@colindices+e2@dimension[2]*(rep(1:e2@dimension[1],diff(e2@rowpointers))-1))))+1
#    z <- .Fortran("_aeexpb_", does not reorder col indices!
#                  as.integer(nrow), as.integer(ncol),
#                  1L,
#                  as.double(e1@entries),  as.integer(e1@colindices),  as.integer(e1@rowpointers),
#                  as.double(e2@entries),  as.integer(e2@colindices),  as.integer(e2@rowpointers),
#                  entries     = vector("double",nzmax),
#                  colindices  = vector("integer",nzmax),
#                  rowpointers = vector("integer",nrow+1),
#                  as.integer(nzmax),
#                  integer(ncol),       double(ncol),
#                  ierr = vector("integer",1),
#                  NAOK=.Spam$NAOK,PACKAGE = "spam"
#                  )
#    if(z$ierr != 0) stop("insufficient space for element-wise exponentiation")
#    nz <- z$rowpointers[nrow+1]-1
#    return(new("spam",entries=z$entries[1:nz],colindices=z$colindices[1:nz],rowpointers=z$rowpointers,
#               dimension=c(nrow,ncol)))
#}

#############################
#getClass("numeric")
#  Extends: "vector"
#getClass("matrix")
#Extends: Class "vector", by class "array", distance 3, with explicit coerce
# Hence we use use vector, especially, to include the NA case that is not of type numeric!
