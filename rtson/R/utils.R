TSON_KIND = "kind"
TSON_SCALAR = -1

add.tson.attribute = function(object, name, value){
  attr = attributes(object)
  if (is.null(attr)){
    attr <- list()
    attr[[name]] <- value
  } else {
    attr[[name]] = value
  }
  attributes(object) <- attr 
  return (object)
}

#' Make a tson map 
#' 
#' Required to generate empty map.
#'
#' @param object A vector or list
#' @return A tson map
#' @export
tson.map = function(object){
  if (is.null(object)) return(NULL)
  return (add.tson.attribute(as.list(object), TSON_KIND, MAP_TYPE))
}

#' Make a tson float32 vector 
#' 
#'
#' @param object A vector or list
#' @return A tson float32 vector 
#' @export
tson.float32.vec = function(object){
  if (is.null(object)) return(NULL)
  return (add.tson.attribute(as.double(object), TSON_KIND, LIST_FLOAT32_TYPE))
}

#' Make a tson int8 vector 
#' 
#'
#' @param object A vector or list
#' @return A tson int8 vector 
#' @export
tson.int8.vec = function(object){
  if (is.null(object)) return(NULL)
  return (add.tson.attribute(as.integer(object), TSON_KIND, LIST_INT8_TYPE))
}

#' Make a tson int16 vector 
#' 
#'
#' @param object A vector or list
#' @return A tson int16 vector 
#' @export
tson.int16.vec = function(object){
  if (is.null(object)) return(NULL)
  return (add.tson.attribute(as.integer(object), TSON_KIND, LIST_INT16_TYPE))
}

#' Make a tson uint8 vector 
#' 
#'
#' @param object A vector or list
#' @return A tson uint8 vector 
#' @export
tson.uint8.vec = function(object){
  if (is.null(object)) return(NULL)
  return (add.tson.attribute(as.integer(object), TSON_KIND, LIST_UINT8_TYPE))
}

#' Make a tson uint16 vector 
#' 
#'
#' @param object A vector or list
#' @return A tson uint16 vector 
#' @export
tson.uint16.vec = function(object){
  if (is.null(object)) return(NULL)
  return (add.tson.attribute(as.integer(object), TSON_KIND, LIST_UINT16_TYPE))
}

#' Make a tson uint32 vector 
#' 
#'
#' @param object A vector or list
#' @return A tson uint32 vector 
#' @export
tson.uint32.vec = function(object){
  if (is.null(object)) return(NULL)
  return (add.tson.attribute(as.integer(object), TSON_KIND, LIST_UINT32_TYPE))
}

#' Make a tson integer
#' 
#'
#' @param object A vector or list
#' @return A tson integer
#' @export
tson.int = function(object){
  if (is.null(object)) return(NULL)
  return (tson.scalar(as.integer(object)))
}

#' Make a tson double
#' 
#'
#' @param object A vector or list
#' @return A tson double
#' @export
tson.double = function(object){
  if (is.null(object)) return(NULL)
  return (tson.scalar(as.double(object)))
}

#' Make a tson character
#' 
#' @param object A vector or list
#' @return A tson character 
#' @export
tson.character = function(object){
  if (is.null(object)) return(NULL)
  return (tson.scalar(as.character(object)))
}

#' Make a tson scalar (ie: singleton)
#' 
#' @param object A vector or list
#' @return A tson scalar 
#' @export
tson.scalar = function(object){
  if (is.null(object)) return(NULL)
  return (add.tson.attribute(object, TSON_KIND, TSON_SCALAR))
}