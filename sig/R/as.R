#' Coerce object to be a siglist
#' 
#' Coerces an object to be a \code{siglist}.
#' @param x Object to coerce.
#' @param ... Passed to other \code{as.siglist} methods.
#' @return An object of class \code{siglist}.
#' @seealso \code{\link{as.sig}}
#' @examples
#' as.siglist(list(
#'   sig(mean),
#'   list(name = "fun", alist(x =,y = 1))
#' ))
#' @export
as.siglist <- function(x, ...)
{
  UseMethod("as.siglist")
}
 
#' @rdname as.siglist
#' @method as.siglist sig
#' @export
as.siglist.sig <- function(x, ...)
{
  sig_list <- list(x)
  names(sig_list) <- x$name
  structure(sig_list, class = c("siglist", "list"))
}
      
#' @rdname as.siglist
#' @method as.siglist list
#' @export
as.siglist.list <- function(x, ...)
{
  structure(
    lapply(x, as.sig),
    class = c("siglist", "list")
  )
}
        
#' @rdname as.siglist
#' @method as.siglist siglist
#' @export
as.siglist.siglist <- function(x, ...)
{
  x
}
 
#' Coerce object to be a sig
#' 
#' Coerces an object to be a \code{sig}.
#' @param x Object to coerce.
#' @param ... Passed to other \code{as.sig} methods.
#' @return An object of class \code{sig}.
#' @seealso \code{\link{as.siglist}}
#' @examples
#' as.sig(
#'   list(name = "fun", alist(x =,y = 1))
#' )
#' @export
as.sig <- function(x, ...)
{
  UseMethod("as.sig")
}
     
#' @rdname as.sig
#' @method as.sig default
#' @export     
as.sig.default <- function(x, ...)
{
  sig(x, ...)
}
          
#' @rdname as.sig
#' @method as.sig siglist
#' @export
as.sig.siglist <- function(x, ...)
{
  if(length(x) > 1)
  {
    warning("Only retaining first sig.")
  }
  x[[1]]
}
        
#' @rdname as.sig
#' @method as.sig list
#' @export
as.sig.list <- function(x, ...)
{
  structure(
    list(
       name = as.character(x$name),
       args = as.list(x$args)
    ),
    class = c("sig", "list")
  )
}
               
#' @rdname as.sig
#' @method as.sig sig
#' @export
as.sig.sig <- function(x, ...)
{
  x
}
 

#' Convert to list
#' 
#' Strips class attributes to return a list.
#' @param x \code{sig}, \code{siglist} or \code{sigreport} object.
#' @param ... Passed from other \code{as.list} methods.
#' @return A list.
#' @examples
#' as.list(sig(read.csv))
#' head(as.list(list_sigs(pkg2env(stats))))
#' as.list(sig_report(baseenv()))
#' @method as.list sig
#' @export
as.list.sig <- function(x, ...)
{
  unclass(x)
}

#' @rdname as.list.sig
#' @method as.list siglist
#' @export
as.list.siglist <- function(x, ...)
{
  unclass(x)
}

#' @rdname as.list.sig
#' @method as.list sigreport
#' @export
as.list.sigreport <- function(x, ...)
{
  unclass(x)
}
