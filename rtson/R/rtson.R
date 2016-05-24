library(R6)
#' Serialize a list 
#'
#' This function convert a list into raw following TSON specification binary-encoded format.
#'
#' @param object A list
#' @return A raw vector
#'
#' @examples
#' ## Example
#' 
#' library(rtson)
#' 
#' list = list(integer=42L,
#'             double=42,
#'             bool=TRUE,
#'             uint8=tson.uint8.vec(c(42,0)),
#'             uint16=tson.uint16.vec(c(42,0)),
#'             uint32=tson.uint32.vec(c(42,0)),
#'             int8=tson.int8.vec(c(42,0)),
#'             int16=tson.int16.vec(c(42,0)),
#'             int32=as.integer(c(42,0)),
#'             float32=tson.float32.vec(c(0.0, 42.0)),
#'             float64=c(42.0,42.0),
#'             map=list(x=42, y=42, label="42"),
#'             list=list("42",42)
#' )
#' 
#' bytes = toTSON(list)
#' @export
toTSON <- function(object) (Serializer$new(object)$bytes)


#' Deserialize a raw vector 
#'
#' This function convert a raw vector into a list following TSON specification binary-encoded format.
#'
#' @param bytes A raw vector
#' @param offset Start offset in raw vector
#' @return A list
#'
#' @examples
#' ## Example
#' 
#' library(rtson)
#' 
#' list = list(integer=42L,
#'             double=42,
#'             bool=TRUE,
#'             uint8=tson.uint8.vec(c(42,0)),
#'             uint16=tson.uint16.vec(c(42,0)),
#'             uint32=tson.uint32.vec(c(42,0)),
#'             int8=tson.int8.vec(c(42,0)),
#'             int16=tson.int16.vec(c(42,0)),
#'             int32=as.integer(c(42,0)),
#'             float32=tson.float32.vec(c(0.0, 42.0)),
#'             float64=c(42.0,42.0),
#'             map=list(x=42, y=42, label="42"),
#'             list=list("42",42)
#' )
#' 
#' bytes = toTSON(list)
#' object = fromTSON(bytes)
#' @export
fromTSON <- function(bytes, offset=NULL) Deserializer$new(bytes, offset)$object


 

 





















