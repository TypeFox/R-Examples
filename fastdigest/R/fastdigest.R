
##' Fast, memory constant hashing of R objects
##'
##' @param obj The object to generate a hash digest for
##' @param ref_serializer (optional) A serializer for reference-style
##' objects, see \code{\link{serialize}}
##'
##' @details \code{obj} will be hashed using R's internal serialization logic
##' with a custom target which applies applying Jenkins' SpookyHash (v2) in a
##' streaming fashion. This avoids (ever) copying the data out of the R object
##' itself, providing both speed and memory constancy.
##'
##' It also guarantees that the "representation" of the R object being hashed
##' is the same as the serialized version would be, if created.
##'
##' @author Gabriel Becker
##' @references Jenkins, B. (2012). SpookyHash: a 128-bit noncryptographic hash.
##' http://burtleburtle.net/bob/hash/spooky.html.
##'
##' @examples
##' x = rnorm(100)
##' h = fastdigest(x)
##' y = x
##' ## force dup
##' y[100] = 1
##' h2 = fastdigest(y)
##' identical(h, h2) #FALSE
##' y[100] = x[100]
##' h3 = fastdigest(y)
##' identical(h, h3) # TRUE
##' @seealso \code{\link{serialize}}
##' @export
fastdigest = function(obj, ref_serializer = NULL)
{
    res = paste(.Call(R_fastdigest, obj, ref_serializer, PACKAGE="fastdigest"),
        collapse="")
    res
}
