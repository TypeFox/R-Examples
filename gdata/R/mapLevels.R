### mapLevels.R
###------------------------------------------------------------------------
### What: Mapping levels
### $Id: mapLevels.R 1991 2015-04-29 03:27:50Z warnes $
### Time-stamp: <2007-04-26 13:16:18 ggorjan>
###------------------------------------------------------------------------

### {{{ mapLevels

###------------------------------------------------------------------------

mapLevels <- function(x, codes=TRUE, sort=TRUE, drop=FALSE,
                      combine=FALSE, ...)
{
  UseMethod("mapLevels")
}

mapLevels.default <- function(x, codes=TRUE, sort=TRUE, drop=FALSE,
                              combine=FALSE, ...)
{
  stop(sprintf("mapLevels can only be used on %s and %s atomic 'x'",
               dQuote("factor"), dQuote("character")))
}

mapLevels.character <- function(x, codes=TRUE, sort=TRUE, drop=FALSE,
                                combine=FALSE, ...)
{
  mapLevels.factor(x=x, codes=codes, sort=sort, drop=drop, ...)
}

## Could coerce character to factor and then use factor method, but that
## is more expensive than simple unique and length used bellow in factor
## method

mapLevels.factor <- function(x, codes=TRUE, sort=TRUE, drop=FALSE,
                             combine=FALSE, ...)
{
  ## --- Argument actions ----

  if(is.factor(x)) { # factor
    if(drop) x <- factor(x)
    nlevs <- nlevels(x)
    levs <- levels(x)
  } else {           # character
    levs <- unique(x)
    nlevs <- length(levs)
    if(sort) levs <- sort(levs, ...)
  }

  ## --- Create a map ---

  map <- vector(mode="list", length=nlevs)
  names(map) <- levs
  if(codes) {
    map[1:nlevs] <- 1:nlevs
  } else {
    map[1:nlevs] <- levs
  }
  class(map) <- "levelsMap"
  map
}

mapLevels.list <- function(x, codes=TRUE, sort=TRUE, drop=FALSE,
                           combine=FALSE, ...)
{
  map <- lapply(x, mapLevels, codes=codes, sort=sort, drop=drop, ...)
  class(map) <- "listLevelsMap"
  if(combine) {
    if(!codes) {
      return(c(map, sort=sort, recursive=TRUE))
    } else {
      stop(sprintf("can not combine integer %s", dQuote("levelsMaps")))
    }
  }
  map
}

mapLevels.data.frame <- function(x, codes=TRUE, sort=TRUE, drop=FALSE,
                                 combine=FALSE, ...)
{
  mapLevels.list(x, codes=codes, sort=sort, drop=drop, combine=combine, ...)
}

### }}}
### {{{ print.*
###------------------------------------------------------------------------

.unlistLevelsMap <- function(x, ind=FALSE)
{
  y <- unlist(x, use.names=FALSE)
  len <- sapply(x, FUN=length)
  names(y) <- rep(names(x), times=len)
  if(ind) {
    return(list(y, rep(1:length(x), times=len), len))
  } else {
    return(y)
  }
}

print.levelsMap <- function(x, ...)
{
  x <- .unlistLevelsMap(x)
  print(x, ...)
}

print.listLevelsMap <- function(x, ...)
{
  class(x) <- "list"
  print(x, ...)
}

### }}}
### {{{ [.*
###------------------------------------------------------------------------

## We need these two since [.list method drops class

"[.levelsMap" <- function(x, i)
{
  classX <- class(x)
  class(x) <- "list"
  x <- x[i]
  class(x) <- classX
  x
}

"[.listLevelsMap" <- function(x, i)
{
  classX <- class(x)
  class(x) <- "list"
  x <- x[i]
  class(x) <- classX
  x
}

### }}}
### {{{ is.*
###------------------------------------------------------------------------

is.levelsMap <- function(x)
  inherits(x=x, what="levelsMap")

is.listLevelsMap <- function(x)
  inherits(x=x, what="listLevelsMap")

.isCharacterMap <- function(x)
{
  if(is(x) == "levelsMap") {
    return(inherits(x=unlist(x), what="character"))
  } else {
    stop(sprintf("can be used only on %s", dQuote("levelsMap")))
  }
}

### }}}
### {{{ as.*
###------------------------------------------------------------------------

as.levelsMap <- function(x, check=TRUE, ...)
{
  if(check)
    .checkLevelsMap(x, method="raw")
  class(x) <- "levelsMap"
  unique(x, ...)
}

as.listLevelsMap <- function(x, check=TRUE)
{
  if(check)
    .checkListLevelsMap(x, method="raw")
  class(x) <- "listLevelsMap"
  x
}

### }}}
### {{{ .check*
###------------------------------------------------------------------------

.checkLevelsMap <- function(x, method) {
  xLab <- deparse(substitute(x))
  also <- "\b"
  if(method == "class") {
    also <- "also"
    if(!is.levelsMap(x))
      stop(sprintf("'%s' must be a %s", xLab, dQuote("levelsMap")))
  }
  if(!is.list(x) || is.null(names(x)))
    stop(sprintf("'%s' must be %s a named list", xLab, also))

  ## Components can be of different length
  ##  if(!all(sapply(x, FUN=length) == 1))
  ##  stop(sprintf("all components of '%s' must have length 1", xLab))
}

.checkListLevelsMap <- function(x, method) {
  xLab <- deparse(substitute(x))
  also <- "\b"
  if(method == "class") {
    also <- "also"
    if(!is.listLevelsMap(x))
      stop(sprintf("'%s' must be a %s", xLab, dQuote("listLevelsMap")))
  }
  if(!is.list(x) || any(!sapply(x, FUN=is.levelsMap)))
    stop(sprintf("'%s' must be %s a list of %s", xLab, also,
                 dQuote("levelsMap")))
  lapply(x, FUN=.checkLevelsMap, method=method)
}

### }}}
### {{{ c.*
###------------------------------------------------------------------------

c.levelsMap <- function(..., sort=TRUE, recursive=FALSE)
{
  x <- list(...)
  class(x) <- "listLevelsMap"
  ## we use recursive=TRUE here because ... is a lists of lists
  c(x, sort=sort, recursive=TRUE)
}

c.listLevelsMap <- function(..., sort=TRUE, recursive=FALSE)
{
  x <- list(...)
  lapply(x, FUN=.checkListLevelsMap, method="class")
  x <- unlist(x, recursive=FALSE)
  if(!recursive) {
    class(x) <- "listLevelsMap"
  } else {
    if(any(!sapply(x, FUN=.isCharacterMap)))
      stop(sprintf("can not combine integer %s", dQuote("levelsMaps")))
    if(!is.null(names(x))) names(x) <- NULL
    x <- unlist(x, recursive=FALSE)
    ## how to merge components with the same name?
    class(x) <- "levelsMap"
    if(sort) x <- sort(x)
    x <- unique(x)
  }
  x
}

### }}}
### {{{ sort
###------------------------------------------------------------------------

sort.levelsMap <- function(x, decreasing=FALSE, na.last=TRUE, ...)
  x[order(names(x), na.last=na.last, decreasing=decreasing)]

### }}}
### {{{ unique
###------------------------------------------------------------------------

unique.levelsMap <- function(x, incomparables=FALSE, ...)
{
  ## Find duplicates
  y <- .unlistLevelsMap(x, ind=TRUE)
  ## Duplicates for values and names combinations
  test <- duplicated(cbind(y[[1]], names(y[[1]])),
                     incomparables=incomparables, ...)
  if(any(test)) {
    if(any(y[[3]] > 1)) { # work with the same structure as in x
      j <- 1
      k <- y[[3]][1]
      empty <- NULL
      for(i in seq(along=x)) { # how slow is this loop?
        tmp <- !test[j:k]
        if(all(!tmp)) { # these components will be empty
          empty <- c(empty, i)
        } else {
          x[[i]] <- x[[i]][tmp]
        }
        j <- j + y[[3]][i]
        k <- k + y[[3]][i + 1]
      }
      if(!is.null(empty))
        x[empty] <- NULL
    } else { # simple one-length components
      x <- x[!test]
    }
  }
  x
}

### }}}
### {{{ mapLevels<-

###------------------------------------------------------------------------

"mapLevels<-" <- function(x, value)
  UseMethod("mapLevels<-")

"mapLevels<-.default" <- function(x, value)
{
  ## --- Checks ---

  classX <- c("integer", "character", "factor")
  if(any(!(class(x) %in% classX)))
    stop(sprintf("'x' must be either: %s", paste(dQuote(classX), collapse=", ")))

  .checkLevelsMap(x=value, method="class")

  ## --- Mapping levels in x ---

  char <- all(sapply(value, is.character))
  int <- all(sapply(value, is.integer))

  if(int) { # codes=TRUE
    if(is.integer(x)) x <- factor(x)
    if(is.factor(x)) levels(x) <- value
    if(is.character(x))
      stop(sprintf("can not apply integer %s to %s",
                   dQuote("levelsMap"), dQuote("character")))
  } else {  # codes=FALSE
    if(!char)
      stop("all components of 'value' must be of the same class")
    if(is.character(x)) x <- factor(x)
    if(is.factor(x)) levels(x) <- value
    if(is.integer(x))
      stop(sprintf("can not apply character %s to %s",
                   dQuote("levelsMap"), dQuote("integer")))
  }
  x
}

"mapLevels<-.list" <- function(x, value)
{
  if(!is.listLevelsMap(value)) {
    if(is.levelsMap(value)) {
      value <- as.listLevelsMap(list(value), check=FALSE)
      ## no need for check as default method does checking anyway
    } else {
      stop(sprintf("'x' must be either %s or %s",
                   dQuote("listLevelsMap"), dQuote("levelsMap")))
    }
  }
  x <- mapply(FUN="mapLevels<-", x=x, value=value, SIMPLIFY=FALSE)
  x
}

"mapLevels<-.data.frame" <- function(x, value)
{
  x[] <- "mapLevels<-.list"(x, value)
  x
}

### }}}
### {{{ Dear Emacs
## Local variables:
## folded-file: t
## End:
### }}}

###------------------------------------------------------------------------
### mapLevels.R ends here
