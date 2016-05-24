
#-----------------------------------------------------------------------------

Ratio <- function(num, den) {
  stopifnot(inherits(num, c("Polyn", "numeric", "integer")))
  stopifnot(inherits(den, c("Polyn", "numeric", "integer")))
  pnum <- as.Polyn(num)
  pden <- as.Polyn(den)
  ratio <- list(numerator=pnum, denominator=pden)
  class(ratio) <- "Ratio"
  ratio
}

#-----------------------------------------------------------------------------

as.Ratio <- function(x, ...) UseMethod("as.Ratio")

as.Ratio.Ratio <- function(x, ...) x

as.Ratio.Polyn <- function(x, ...) Ratio(x, 1)

as.Ratio.default <- function(x, ...) as.Ratio(as.Polyn(x, ...), ...)

#-----------------------------------------------------------------------------

Rnumerator <- function(r) {
  stopifnot(inherits(r, "Ratio"))
  r$numerator
}
Rdenominator <- function(r) {
  stopifnot(inherits(r, "Ratio"))
  r$denominator
}

#-----------------------------------------------------------------------------

as.character.Ratio <- function(x, ...) {
  paste("(", as.character(Rnumerator(x), ...), ") / (", 
    as.character(Rdenominator(x), ...), ")", sep="")
}

print.Ratio <- function(x, ...) {
  print(as.character(x, ...))
  invisible(x)
}

#-----------------------------------------------------------------------------
# Arithmetic operators

`+.Ratio` <- function(r1, r2) {
  stopifnot(inherits(r1, c("Ratio", "Polyn", "numeric", "integer")))
  stopifnot(inherits(r2, c("Ratio", "Polyn", "numeric", "integer")))
  r1<-as.Ratio(r1)
  r2<-as.Ratio(r2)
  Ratio(r1$numerator*r2$denominator+r2$numerator*r1$denominator, 
    r1$denominator*r2$denominator)
}

`-.Ratio` <- function(r1, r2) {
  stopifnot(inherits(r1, c("Ratio", "Polyn", "numeric", "integer")))
  stopifnot(inherits(r2, c("Ratio", "Polyn", "numeric", "integer")))
  r1<-as.Ratio(r1)
  r2<-as.Ratio(r2)
  Ratio(r1$numerator*r2$denominator-r2$numerator*r1$denominator, 
    r1$denominator*r2$denominator)
}

`*.Ratio` <- function(r1, r2) {
  stopifnot(inherits(r1, c("Ratio", "Polyn", "numeric", "integer")))
  stopifnot(inherits(r2, c("Ratio", "Polyn", "numeric", "integer")))
  r1<-as.Ratio(r1)
  r2<-as.Ratio(r2)
  Ratio(r1$numerator*r2$numerator, r1$denominator*r2$denominator)
}

`/.Ratio` <- function(r1, r2) {
  stopifnot(inherits(r1, c("Ratio", "Polyn", "numeric", "integer")))
  stopifnot(inherits(r2, c("Ratio", "Polyn", "numeric", "integer")))
  r1<-as.Ratio(r1)
  r2<-as.Ratio(r2)
  Ratio(r1$numerator*r2$denominator, r1$denominator*r2$numerator)
}

`^.Ratio` <- function(r, n) {
  stopifnot(inherits(r, "Ratio"))
  stopifnot(inherits(n, c("numeric", "integer")))
  stopifnot(length(n)==1)
  stopifnot(round(n)==n & n>=0)
  if(n==0) as.Ratio(as.Polyn(1, 0))
  else if (n==1) r
  else Ratio(r$numerator^n, r$denominator^n)
}

#-----------------------------------------------------------------------------
# Relational operators

`==.Ratio` <- function(r1, r2) identical(r1, r2)

`!=.Ratio` <- function(r1, r2) !identical(r1, r2)

`<.Ratio` <- function(r1, r2) as.logical(NA)
`<=.Ratio` <- function(r1, r2) as.logical(NA)
`>.Ratio` <- function(r1, r2) as.logical(NA)
`>=.Ratio` <- function(r1, r2) as.logical(NA)

#-----------------------------------------------------------------------------
