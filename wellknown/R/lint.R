#' Validate WKT strings
#'
#' @export
#' @param str A WKT string
#' @return A logical (\code{TRUE} or \code{FALSE})
#' @examples
#' lint("POINT (1 2)")
#' lint("POINT (1 2 3)")
#' lint("LINESTRING EMPTY")
#' lint("LINESTRING (100 0, 101 1)")
#' lint("MULTIPOINT EMPTY")
#' lint("MULTIPOINT ((1 2), (3 4))")
#' lint("MULTIPOINT ((1 2), (3 4), (-10 100))")
#' lint("POLYGON ((1 2, 3 4, 0 5, 1 2))")
#' lint("POLYGON((20.3 28.6, 20.3 19.6, 8.5 19.6, 8.5 28.6, 20.3 28.6))")
#' lint("MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, 40 10, 10 20, 5 10, 15 5)))")
#' lint("TRIANGLE ((0 0, 0 1, 1 1, 0 0))")
#' lint("TRIANGLE ((0.1 0.1, 0.1 1.1, 1.1 1.1, 0.1 0.1))")
#' lint("CIRCULARSTRING (1 5, 6 2, 7 3)")
#' lint("CIRCULARSTRING (1 5, 6 2, 7 3, 5 6, 4 3)")
#' lint('COMPOUNDCURVE (CIRCULARSTRING (1 0, 0 1, -1 0), (-1 0, 2 0))')
lint <- function(str) {
  type <- get_type(str, ignore_case = TRUE)
  if (length(type) == 0) return(FALSE)
  str <- str_trim_(gsub(toupper(type), "", str))
  rule <- switch(type,
                 Point = rule_point,
                 Linestring = rule_linestring,
                 Multipoint = rule_multipoint,
                 Polygon = rule_polygon,
                 Multipolygon = rule_multipolygon,
                 Triangle = rule_triangle,
                 Circularstring = rule_circularstring,
                 Compoundcurve = rule_compoundcurve)
  all(unlist(Map(function(x, y) eval(parse(text = x))(str, y), names(rule), rule)))
}

# collapse rules into one string
coll <- function(x) paste0(x, collapse = "")

# grepl with concatenating rule strings together first
grepl2 <- function(str, x) {
  grepl(coll(x), str)
}

# put rules in vector separated by OR's
vor <- function(...) {
  rlz <- list(...)
  unlist(Map(function(x,y) c(x,y), rlz, c(rep("|", length(rlz) - 1), "") ))
}

# use grepl on many tests, and give logical for all tests
greplall <- function(rlz, str) {
  all(sapply(rlz, function(z) grepl(coll(z), x = str)))
}

# make a string a repeat pattern
repeat_ <- function(x) sprintf("(%s)*", coll(x))

# make a string repeat a pattern n or more times
repeat_n <- function(x, n = 3) sprintf("(%s){%s,}", coll(x), n)

# short-hand nouns
number <- "[+-]?(\\d*\\.)?\\d+"
space <- "\\s"
spaceif <- "\\s?"
comma <- ","
lp <- "^\\("
lp_ <- "\\("
rp <- "\\)$"
rp_ <- "\\)"
empty <- "^EMPTY$"
rule_point_empty <- rule_multipoint_empty <- rule_linestring_empty <- rule_polygon_empty <- rule_multipolygon_empty <- empty
pt <- c(number, space, number)
pt3 <- c(number, space, number, space, number)
commapt <- c(comma, spaceif, pt)
multipt <- c(lp_, pt, rp_)
linestr <- c(lp_, pt, repeat_(commapt), rp_)
polygonstr <- c(lp_, linestr, rp_)
commapolygon <- c(comma, spaceif, polygonstr)
reppolygonstr <- c(lp_, linestr, rp_, repeat_(commapolygon))
commamultipt <- c(comma, spaceif, multipt)
commalinestr <- c(comma, spaceif, linestr)
commapt3 <- c(comma, spaceif, pt3)
trep_point <- 'POINT'
trep_polygon <- 'POLYGON'
trep_multipolygon <- 'MULTIPOLYGON'
trep_linestring <- 'LINESTRING'
trep_multilinestring <- 'MULTILINESTRING'
trep_triangle <- 'TRIANGLE'
trep_circularstring <- 'CIRCULARSTRING'

# point rules
rule_point_2d <- c(lp, pt, rp)
rule_point_3d <- c(lp, number, space, number, space, number, rp)
rule_point_4d <- c(lp, number, space, number, space, number, space, number, rp)
rule_point <- list(grepl2 = vor(empty, rule_point_2d, rule_point_3d, rule_point_4d))

# multipoint rules
rule_multipoint_2d <- c(lp, multipt, repeat_(commamultipt), rp)
rule_multipoint <- list(grepl2 = vor(empty, rule_multipoint_2d))

# linestring rules
rule_linestring_2d <- c(lp, pt, repeat_(commapt), rp)
rule_linestring_3d <- c(lp, pt3, repeat_(commapt3), rp)
rule_linestring <- list(grepl2 = vor(empty, rule_linestring_2d, rule_linestring_3d))

# multilinestring rules
# xxxxx

# polygon rules
rule_polygon_2d <- c(lp, linestr, repeat_(commalinestr), rp)
rule_polygon <- list(grepl2 = vor(empty, rule_polygon_2d))

# multipolygon rules
rule_multipolygon_2d <- c(lp, reppolygonstr, rp)
rule_multipolygon <- list(grepl2 = vor(empty, rule_multipolygon_2d))

# triangle rules
rule_triangle_2d <- c(lp, lp_, pt, repeat_(commapt), rp_, rp)
rule_triangle <- list(grepl2 = vor(empty, rule_triangle_2d))

# circularstring rules
rule_circularstring_2d <- c(lp, pt, repeat_(commapt), rp)
rule_circularstring <- list(grepl2 = vor(empty, rule_circularstring_2d))
# rule_circularstring <- list(grepl2 = vor(empty, rule_circularstring_2d), check_circularstring = NULL)

# compoundcurve rules
rule_compoundcurve_2d <- c(lp, trep_circularstring,  lp_, pt, repeat_(commapt), rp_, rp)
rule_compoundcurve <- list(grepl2 = vor(empty, rule_compoundcurve_2d))

check_circularstring <- function(x, ...) {
  length(strsplit(x, ",")[[1]]) == 3
}

# 'COMPOUNDCURVE EMPTY'
# 'COMPOUNDCURVE (CIRCULARSTRING (1 0, 0 1, -1 0), (-1 0, 2 0))'
# 'COMPOUNDCURVE (CIRCULARSTRING (1 1, 1 1, 1 1), (1 1, 3 5, 5 4))'
