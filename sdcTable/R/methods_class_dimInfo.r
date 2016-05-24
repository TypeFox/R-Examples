########################################
### methods only for class 'dimInfo' ###
########################################
#' @aliases get.dimInfo,dimInfo,character-method
#' @rdname get.dimInfo-method
setMethod(f="get.dimInfo", signature=c("dimInfo", "character"),
  definition=function(object, type) {
    if ( !type %in% c("strInfo", "dimInfo", "varName", "strID", "posIndex") ) {
      stop("get.dimInfo:: argument 'type' is not valid!\n")
    }
    if ( type == "strInfo" ) {
      return(g_str_info(object))
    }
    if ( type == "dimInfo" ) {
      return(g_dim_info(object))
    }
    if ( type == "varName" ) {
      return(g_varname(object))
    }
    if ( type == "strID" ) {
      return(g_str_id(object))
    }
    if ( type == "posIndex" ) {
      return(g_pos_index(object))
    }
  }
)

#' @aliases set.dimInfo,dimInfo,character,character-method
#' @rdname set.dimInfo-method
setMethod(f="set.dimInfo", signature=c("dimInfo", "character", "character"),
  definition=function(object, type, input) {
    if ( !type %in% c("strID") ) {
      stop("set.dimInfo:: check argument 'type'!\n")
    }
    if ( type == "strID" ) {
      s_str_id(object) <- input
    }
    validObject(object)
    return(object)
  }
)

setMethod(f="g_str_info", signature=c("dimInfo"), definition=function(object) {
  return(object@strInfo)
})
setMethod(f="g_dim_info", signature=c("dimInfo"), definition=function(object) {
  return(object@dimInfo)
})
setMethod(f="g_varname", signature=c("dimInfo"), definition=function(object) {
  return(object@vNames)
})
setMethod(f="g_str_id", signature=c("dimInfo"), definition=function(object) {
  return(object@strID)
})
setMethod(f="g_pos_index", signature=c("dimInfo"), definition=function(object) {
  return(object@posIndex)
})
setReplaceMethod(f="s_str_id", signature=c("dimInfo", "character"), definition=function(object, value) {
  object@strID <- value
  return(object)
})

