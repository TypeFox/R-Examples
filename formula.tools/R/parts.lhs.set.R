# -----------------------------------------------------------------------------
# REPLACEMENT : lhs<-
# -----------------------------------------------------------------------------

#' @rdname formula.parts
#' @name lhs<-
#' @export lhs<-
#' @include parts.lhs.get.R

setGeneric( 'lhs<-', function( x, value ) standardGeneric('lhs<-') )


# -------------------------------------
# SINGLULAR: call, formula
# -------------------------------------

#' @rdname formula.parts
#' @aliases .replace.lhs.singular

.replace.lhs.singular <-  function( x, value ) {
  x[[2]] <- value 
  x 
}


#' @rdname formula.parts
#' @name lhs<-
#' @aliases lhs<-,call-method
setReplaceMethod( 'lhs', 'call', .replace.lhs.singular )


#' @name lhs<-
#' @rdname formula.parts
#' @aliases lhs<-,formula-method
setReplaceMethod( 'lhs', 'formula' , .replace.lhs.singular )




# -------------------------------------
# LIST AND VECTORS: expression, list
# -------------------------------------
# .replace.lhs.plural <- function( x, value ) {
# 
#     if( length(value) == 1 ) {
#       for( i in 1:length(x) ) lhs( x[[i]] ) <- value 
#     } else {  
#       if( length(x) != length(value) ) 
#         stop( "Cannot change the lhs.  Arguments have different lengths" )
# 
#       for( i in 1:length(x) ) lhs(x[[i]] ) <- value[[i]]
#     }
# 
#     x
# }        


#' @rdname formula.parts
#' @aliases .replace.lhs.plural

.replace.lhs.plural <- function( x, value ) { 
  
  if( length(value) == 1 ) { 
    for( i in 1:length(x) ) lhs( x[[i]] ) <- value 
    
  } else if( length(x) == length(value) ) {
    for( i in 1:length(x) ) lhs( x[[i]] ) <- value[[i]]
    
  } else { 
    warning( "length of object != length of lhs replacement" )
  }
  
  x 
  
}

#' @name lhs<-  
#' @rdname formula.parts 
#' @aliases lhs<-,expression-method
setReplaceMethod( 'lhs', c('expression','ANY') , .replace.lhs.plural )

#' @name lhs<-
#' @rdname formula.parts
#' @aliases lhs<-,list-method 
setReplaceMethod( 'lhs', 'list' , .replace.lhs.plural )

