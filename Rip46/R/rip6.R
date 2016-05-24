#' IPv6 Conversion Functions
#' 
#'  Functions to manipulate objects of class \code{"ip6"} representing IP version 6 addresses.
#'
#' @name ip6
#' @param x An object to be converted.
NULL 

#' @export
#' @rdname ip6
as.ip6 <- function(x) UseMethod("as.ip6")

#' @export
as.ip6.ip6 <- identity

#' @export
#' @rdname ip6
#' @examples as.ip6("DE:AD:BE:EF:01:23:45:67")
as.ip6.character <- hostToIp6

#' @export
#' @rdname ip6
#' @param ...   further arguments passed to or from other methods.
as.character.ip6 <- function(x, ...) ip6ToHost(x)

#' @export
as.data.frame.ip6 <- as.data.frame.complex



#' @export
format.ip6 <- function(x, ...) ip6ToHost(x)

#' @export
print.ip6 <- function(x, ...) print(ip6ToHost(x), ...)

#' @export
`[.ip6` <- function(x, ...) `oldClass<-`( unclass(x)[...], oldClass(x))

#' @export
`[<-.ip6` <- function(x,i, value) {
  `oldClass<-`( `[<-`(unclass(x), i, as.ip6(value)), oldClass(x))
}

# Make ip4 work with colClasses
setClass("ip6")
setAs("character","ip6", function(from) hostToIp4(from) )

#' @export
Ops.ip6 <- function(e1, e2) 
  stop(gettextf("%s not defined for \"ip6\" objects", .Generic), domain = NA)

#' @export
Math.ip6 <- function(x, ...) 
  stop(gettextf("%s not defined for \"ip6\" objects", .Generic), domain = NA)