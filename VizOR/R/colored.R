##' Create a colored factor,
##' an object of class \code{c('colored','factor')} or \code{c('colored','ordered')},
##' from a character vector or factor \code{x}.
##' 
##' Colored factors permit persistent association of color to a factor,
##' supporting consistent graphical treatment across multiple plots.
##' 
##' @param x A character vector or factor
##' @param color.key A named vector or list mapping factor levels (the
##' names) to colors (the values). If missing (as is likely in initial
##' exploratory analyses), this is generated automatically as a
##' convenience to the user. In case an ordered factor is generated,
##' the order of the levels is determined by the ordering in
##' \code{color.key}.
##' @param ordered A logical value stating whether the factor should be ordered
##' @param default An optional argument; when provided, it gives a string in
##'   names(color.key) to which all non-matching values of x should be mapped.
##'   The effect will be to collect possibly many factor levels under a single
##'   key entry.  A common choice might be, for example, \code{default='OTHER'}.
##' @return An object of class \code{c('colored',['ordered',]'factor')}
##' @author David C. Norris
##' @keywords classes category
##' @examples
##' # This is an example of de novo construction of a colored factor
##' weekdays <- colored(c("Mon","Tue","Wed","Thu","Fri"),
##'                     color.key=c(Mon="blue",Tue="red",Wed="yellow",
##'                                 Thu="purple", Fri="green"))
##' # This demonstrates how one might use the 'colored' constructor
##' # to expand the level set of an existing factor.
##' week <- colored(weekdays,
##'                 color.key=c(Sun="white", key(weekdays), Sat="gray"),
##'                 ordered=TRUE)
##'
##' # Note that 'droplevels.factor' works fine on colored factors
##' levels(week)
##' levels(droplevels(week))
##' @export colored
colored <- function(x, color.key, ordered=is.ordered(x), default=NA){
  ordered <- eval(ordered) # don't get fooled by lazy evaluation of this default!
  if(missing(color.key)){
    if(is(x,"colored")){
      color.key <- key(x)[names(key(x)) %in% unique(x)] # so ordered(.) drops unused levels, as factor(.) does
    } else {
      warning("Generating colored factor with arbitrary colors, as 'color.key' unspecified")
      levels <- if(is.factor(x)) levels(x) else as.character(unique(x))
      color.key <- rainbow(length(levels))
      names(color.key) <- levels
    }
  }
  ## TODO: If missing(default), assert all(x %in% names(color.key))?
  x <- match(x, names(color.key), nomatch=match(default, names(color.key)))
  x <- factor(names(color.key)[x], levels=names(color.key), ordered=ordered)
  attr(x, 'colors') <- unname(unlist(color.key))
  class(x) <- c('colored', class(x))
  x
}
