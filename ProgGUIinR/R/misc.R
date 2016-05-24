## Miscellaneous functions

faux_require <- function(pkg) {
  do.call("require", list(package=pkg))
}

##' List available data frames
##'
##' @param envir the environment to search in
##' @return a character vector with data frame names
##' @export
avail_dfs <- function(envir=parent.frame()) {
  x <- ls(envir=envir)
  x[sapply(mget(x, envir=envir), is.data.frame)]
}



##' return variable names in a data frame or environment matching some condition
##'
##' @param x a data frame,data frame name in the global workspace or an environment
##' @param cond a function to apply to each column of the data frame Only those matching a returned
##' @return a character vector of variable names.
##' @export
##' @examples
##' library(MASS) # for Cars93 data frame
##' find_vars(Cars93, cond=is.factor)
find_vars <- function(x, cond=function(x) TRUE) UseMethod("find_vars")

##' default method
##'
##' @S3method find_vars default
##' @method find_vars default
##' @rdname find_vars
find_vars.default <- function(x, cond=function(x) TRUE) character(0)

##' character method
##' @S3method find_vars character
##' @method find_vars character
##' @rdname find_vars
find_vars.character <- function(x, cond=function(x) TRUE) {
  ## df is data frame or name of data frame in given environment
  x <- get(x, envir=.GlobalEnv)
  find_vars(x, cond)
}

##' data frame method
##'
##' @S3method find_vars data.frame
##' @method find_vars data.frame
##' @rdname find_vars
find_vars.data.frame <- function(x, cond=function(x) TRUE) {
  names(x)[sapply(x, cond)]
}

##' find variables
##'
##' @S3method find_vars environment
##' @method find_vars environment
##' @rdname find_vars
find_vars.environment <- function(x, cond=function(x) TRUE) {
  ind <- sapply(ls(x), function(i) cond(get(i, envir=x)))
  ls(x)[ind]
}


##' sprintf with template depending on integer valued n
##'
##' @param n number passed to sprintf call
##' @param msg1 \code{\link{sprintf}} template when \code{n == 1}. 
##' @param msg2 \code{\link{sprintf}} template when \code{n > 1}.
##' @param ... passed to \code{ngettext}
##' @return a character string
##' @export
nsprintf <- function(n, msg1, msg2,...)
  ngettext(n, sprintf(msg1, n), sprintf(msg2, n), ...)

##################################################
## Date functions

##' is this a valid date?
##'
##' @param year year as number
##' @param month Month in 1 to 12
##' @param day in 1:12
##' @return logical
##' @export
validDate <- function(year, month, day) 
  !is.na(as.Date(sprintf("%s-%s-%s", year, month, day), format="%Y-%m-%d"))

##' Day of week for a given year, month day
##'
##' from chron with slight change to arguments
##' @param year year
##' @param month month as integer 1 to 12
##' @param day day of month
##' @return day of week in 0:6
##' @export
day.of.week <- function (year, month, day) {
    ix <- year + trunc((month - 14)/12)
    jx <- (trunc((13 * (month + 10 - (month + 10)%/%13 * 12) - 
        1)/5) + day + 77 + (5 * (ix - (ix%/%100) * 100))%/%4 + 
        ix%/%400 - (ix%/%100) * 2)
    jx%%7
}

##' how many days in a month
##'
##' @param year year as integer
##' @param month month in 1 to 12
##' @return number of days in the given month
##' @export
days.in.month <- function(year, month) {
  for(i in c(31, 30, 29, 28)) {
    if(validDate(year, month, i))
      return(i)
  }
}

##' 0-based week of month
##'
##' @param year year
##' @param month month in 1 .. 12
##' @param day day of year
##' @return week (0-based) of the given day
##' @export
week.of.month <- function(year, month, day) {
  first.day <- day.of.week(year, month, 1)
  (first.day + day - 1) %/% 7
}



##' icon functions
##'
##' Idea is taken from ggplot2
##' @param x object
##' @param ... dots
##' @export
##' @return a recipe to draw an icon
##' @examples
##' x <- rnorm(10)
##' library(grid); grid.newpage(); grid.draw(make_icon(x))
make_icon <- function(x, ...) UseMethod("make_icon")

##' method for numeric class
##'
##' @inheritParams make_icon
##' @S3method make_icon numeric
##' @method make_icon numeric
##' @rdname make_icon
make_icon.numeric <- function(x, ...) grid::textGrob("#", gp=grid::gpar(cex=1.2, col="gray50"))

##' method for numeric class
##'
##' @inheritParams make_icon
##' @S3method make_icon factor
##' @method make_icon factor
##' @rdname make_icon
make_icon.factor <- function(x, ...)  grid::textGrob("f", gp=grid::gpar(cex=1.2, col="gray50", fontface="italic"))

##' method for numeric class
##'
##' @inheritParams make_icon
##' @S3method make_icon character
##' @method make_icon character
##' @rdname make_icon
make_icon.character <- function(x, ...) {
  cex <- 1.2
  grid::grid.draw(grid::textGrob("a", x=.25, y=.25, gp=grid::gpar(cex=cex, col="gray70")))
  grid::grid.draw(grid::textGrob("b", x=.5, y=.75, gp=grid::gpar(cex=cex, col="gray60")))
  grid::grid.draw(grid::textGrob("c", x=.75, y=.5, gp=grid::gpar(cex=cex, col="gray50")))
}

##' method for numeric class
##'
##' @inheritParams make_icon
##' @S3method make_icon logical
##' @method make_icon logical
##' @rdname make_icon
make_icon.logical <- function(x, ...) {
  cex <- 1.2
  grid::grid.draw(grid::textGrob("FALSE", x=.75, y=.25, gp=grid::gpar(cex=cex, col="gray50")))
  grid::grid.draw(grid::textGrob("TRUE", x=.25, y=.75, gp=grid::gpar(cex=cex, col="gray70")))
}

##' method for numeric class
##'
##' @inheritParams make_icon
##' @S3method make_icon data.frame
##' @method make_icon data.frame
##' @rdname make_icon
make_icon.data.frame <- function(x, ...) {
  cex <- 1.2
  grid::grid.draw(grid::textGrob("data.frame", gp=grid::gpar(cex=cex, col="gray70")))
}

##' make an icon for type of graph
##'
##' Stolen from ggplot -- but this is now deprecated!
##' @param x type of graph. E.g., "boxplot", "abline", ...
##' @param ... dots
##' @return a recipe to draw an icon
##' @export
##' @examples
##' ##' ## Doesn't work as of ggplot2 version 0.9.2
##' ##library(grid); grid.newpage(); grid.draw(make_graph_icon("hist"))
##'
make_graph_icon <- function(x,...) {
    .Deprecated("NA")
##   require(ggplot2)
##   switch(x,
##          "abline"=ggplot2::GeomAbline$icon(),
##          "boxplot"=ggplot2::GeomBoxplot$icon(),
##          "bwplot"=ggplot2:::GeomBoxplot$icon(),
##          "contour"=ggplot2:::GeomContour$icon(),
##          "density"=ggplot2:::GeomDensity$icon(),
##          "hist"=ggplot2:::GeomHistogram$icon(),
##          "histogram"=ggplot2:::GeomHistogram$icon(),
##          "lines"=ggplot2:::GeomLine$icon(),
##          "points"=ggplot2:::GeomPoint$icon(),
##          "polygon"=ggplot2:::GeomPolygon$icon(),
##          "qqplot"=ggplot2:::GeomQuantile$icon(),
##          "qqnorm"=ggplot2:::GeomQuantile$icon(),
##          "rect" = ggplot2:::GeomRect$icon(),
##          "rug" = ggplot2:::GeomRug$icon(),
##          ggplot2:::GeomPoint$icon()
##        )
}



##' Lorem Ipsum text
##'
##' The Lorem ipsum text is used in typesetting to fill space when showing a mockup.
##' @format chr Lorem ipsum dolor sit amet, consectetur adipisicing elit, ...
##' @export
LoremIpsum <- "
Lorem ipsum dolor sit amet, consectetur adipisicing elit,
sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi
ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit
in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur
sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit
anim id est laborum.

Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque
laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi
architecto beatae vitae dicta sunt explicabo. Nemo enim ipsam voluptatem quia voluptas
sit aspernatur aut odit aut fugit, sed quia consequuntur magni dolores eos qui ratione
voluptatem sequi nesciunt. Neque porro quisquam est, qui dolorem ipsum quia dolor sit amet,
consectetur, adipisci velit, sed quia non numquam eius modi tempora incidunt ut labore et
dolore magnam aliquam quaerat voluptatem. Ut enim ad minima veniam, quis nostrum
exercitationem ullam corporis suscipit laboriosam, nisi ut aliquid ex ea commodi
consequatur? Quis autem vel eum iure reprehenderit qui in ea voluptate velit esse quam
nihil molestiae consequatur, vel illum qui dolorem eum fugiat quo voluptas nulla pariatur?
"
