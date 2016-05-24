# -----------------------------------------------------------------------------
# rhs
#   extract and manipulate the right-hand side of R objects
# -----------------------------------------------------------------------------

#' @include parts.R
NULL

# -----------------------------------------------------------------------------
# REPLACEMENT rhs<-
# -----------------------------------------------------------------------------


#' @name rhs<-
#' @aliases rhs<-
#' @rdname formula.parts
#' @export rhs<-

  setGeneric( 'rhs<-', function(x,value) standardGeneric('rhs<-') )

# -------------------------------------
# SINGULAR: call, formula
# -------------------------------------

#' @rdname formula.parts
#' @aliases .replace.rhs.singular
.replace.rhs.singular <-  function( x, value ) {
  x[[3]] <- value 
  x 
}                                                    

#' @rdname formula.parts
#' @name rhs<- 
#' @aliases rhs<-,call-method
setReplaceMethod( 'rhs', 'call' , .replace.rhs.singular )

#' @rdname formula.parts
#' @name rhs<- 
#' @aliases rhs<-,formula-method
setReplaceMethod( 'rhs', 'formula' , .replace.rhs.singular )



# -------------------------------------
# PLURAL: LIST AND VECTORS: expression, list
# 
#  Note: 
#   - It is possible to have the rhs contain more than
#     one value, e.g. rhs(e) <- 1:3.  Because of the 
#     ambiguity, we do not do multiple replaces.
# -------------------------------------
# .replace.rhs.plural <- function( x, value ) {
# 
#     if( length(value) == 1 ) {
#       for( i in 1:length(x) ) rhs( x[[i]] ) <- value 
# 
#     } else {  
# 
#       if( length(x) != length(value) ) 
#         stop( "Cannot change the rhs. Arguments have different lengths." )
# 
#       for( i in 1:length(x) ) rhs( x[[i]] ) <- value[[i]]
# 
#     }
# 
#     x
# }        

#' @rdname formula.parts
#' @aliases .replace.ths.plural
.replace.rhs.plural <- function( x, value ) {

  if( length(value) == 1 ) { 
    for( i in 1:length(x) ) rhs( x[[i]] ) <- value 
    
  } else if( length(x) == length(value) ) {
    for( i in 1:length(x) ) rhs( x[[i]] ) <- value[[i]]
    
  } else { 
    warning( "length of object != length of rhs replacement" )
  }
  
  x     
  
}

#' @name rhs<-
#' @rdname formula.parts
#' @aliases rhs<-,expression-method
setReplaceMethod( 'rhs', 'expression' , .replace.rhs.plural )

#' @name rhs<-
#' @rdname formula.parts
#' @aliases rhs<-,list-method
setReplaceMethod( 'rhs', 'list' , .replace.rhs.plural )

