##' Reorder levels of a colored factor.
##' 
##' @aliases relevel
##' @param x A colored factor
##' @param ref The reference level
##' @param \dots Unused
##' @return A colored factor
##' @method relevel colored
##' @S3method relevel colored
##' @seealso \code{\link{relevel}}
##' @author David C. Norris
##' @keywords category color
##' @examples # TODO: Provide an example
##' @export relevel
relevel.colored <- function(x, ref, ...){
  colored(NextMethod("relevel"),
          color.key=key(x)[levels(relevel(factor(NA, levels=levels(x)), ref))])
}
