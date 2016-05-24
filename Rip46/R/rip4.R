#' IPv4 Conversion Functions
#' 
#'  Functions to manipulate objects of class \code{"ip4"} representing IP version 4 addresses.
#'
#' @name ip4
#' @param x An object to be converted.
#' @import methods
NULL 

#' @rdname ip4
#' @export
as.ip4 <- function(x) UseMethod("as.ip4")

#' @rdname ip4
#' @export
as.ip4.numeric <- mySqlToIp4

#' @export
as.ip4.ip4 <- identity

#' @rdname ip4
#' @export
#' @examples as.ip4("192.168.0.1")
#' @examples read.table(text="128.0.0.1", colClasses = "ip4")
as.ip4.character <- hostToIp4

#' @export
as.data.frame.ip4 <- as.data.frame.numeric

#' @rdname ip4
#' @param ...   further arguments passed to or from other methods.
#' @export
as.character.ip4 <- function(x, ...) ip4ToHost(x)

#' @export
format.ip4 <- function(x, ...) ip4ToHost(x)

#' @export
print.ip4 <- function(x, ...) print(ip4ToHost(x), ...)

#' @export
`[.ip4` <- function(x, ...) `oldClass<-`( unclass(x)[...], oldClass(x))

#' @export
`[<-.ip4` <- function(x,i, value) {
  `oldClass<-`( `[<-`(unclass(x), i, as.ip4(value)), oldClass(x))
}

# Make ip4 work with colClasses
setClass("ip4")
setAs("character","ip4", function(from) hostToIp4(from) )

#' @export
Ops.ip4 <- function(e1, e2) 
  stop(gettextf("%s not defined for \"ip4\" objects", .Generic), domain = NA)

#' @export
Math.ip4 <- function(x, ...) 
  stop(gettextf("%s not defined for \"ip4\" objects", .Generic), domain = NA)