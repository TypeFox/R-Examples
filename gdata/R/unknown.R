### unknown.R
###------------------------------------------------------------------------
### What: Change given unknown value to NA and vice versa
### $Id: unknown.R 1797 2014-04-05 18:19:49Z warnes $
### Time-stamp: <2007-04-26 13:16:10 ggorjan>
###------------------------------------------------------------------------

### {{{ isUnknown

###------------------------------------------------------------------------

isUnknown <- function(x, unknown=NA, ...)
  UseMethod("isUnknown")

isUnknown.default <- function(x, unknown=NA, ...)
{
  if(is.list(unknown)) unknown <- unlist(unknown)
  ret <- x %in% unknown
  if(any(is.na(unknown))) ret <- ret | is.na(x)
  ret
}

isUnknown.POSIXlt <- function(x, unknown=NA, ...)
{
  ## FIXME: codetools say
  ## isUnknown.POSIXlt: wrong number of arguments to as.character
  if(is.list(unknown) && !inherits(x=unknown, what="POSIXlt")) {
    unknown <- lapply(unknown, FUN=as.character, ...)
  } else {
    unknown <- as.character(x=unknown, ...)
  }

  if(is.list(x) && !inherits(x=x, what="POSIXlt")) {
    x <- lapply(x, FUN=as.character, ...)
  } else {
    x <- as.character(x=x, ...)
  }

  isUnknown.default(x=as.character(x), unknown=as.character(unknown))
}

isUnknown.list <- function(x, unknown=NA, ...) {
  unknown <- .unknownList(x=x, unknown=unknown)
  x <- mapply(FUN="isUnknown", x=x, unknown=unknown, ..., SIMPLIFY=FALSE)
  x
}

isUnknown.data.frame <- function(x, unknown=NA, ...)
{
  x[] <- isUnknown.list(x, unknown=unknown, ...)
  x
}

isUnknown.matrix <- function(x, unknown=NA, ...)
  apply(X=x, MARGIN=ifelse(ncol(x) > nrow(x), 1, 2), FUN=isUnknown,
        unknown=unknown)

### }}}
### {{{ unknownToNA

###------------------------------------------------------------------------

unknownToNA <- function(x, unknown, warning=FALSE, ...)
  UseMethod("unknownToNA")

unknownToNA.default <- function(x, unknown, warning=FALSE, ...)
{
  if(warning) {
    if(any(is.na(x)))
      warning("'x' already has NA")
  }
  is.na(x) <- isUnknown(x=x, unknown=unknown)
  x
}

unknownToNA.factor <- function(x, unknown, warning=FALSE, ...)
{
  ## could put this func into default method, but I need unlisted unknown
  ## for levels handling
  if(warning) {
    if(any(is.na(x)))
      warning("'x' already has NA")
  }
  if(is.list(unknown)) unknown <- unlist(unknown)
  ## Levels handling - read help page on this
  levs <- levels(x)
  levs <- levs[!(levs %in% unknown)]
  factor(x, levels=levs)
}

unknownToNA.list <- function(x, unknown, warning=FALSE, ...)
{
  unknown <- .unknownList(x=x, unknown=unknown)
  x <- mapply(FUN="unknownToNA", x=x, unknown=unknown, warning=warning,
              SIMPLIFY=FALSE)
  return(x)
}

unknownToNA.data.frame <- function(x, unknown, warning=FALSE, ...)
{
  x[] <- unknownToNA.list(x=x, unknown=unknown, warning=warning)
  x
}

### }}}
### {{{ NAToUnknown

###------------------------------------------------------------------------

NAToUnknown <- function(x, unknown, force=FALSE, call.=FALSE, ...)
  UseMethod("NAToUnknown")

NAToUnknown.default <- function(x, unknown, force=FALSE, call.=FALSE, ...)
{
  if(length(as.character(unknown)) != 1) # as.character allows also POSIXlt
    stop("'unknown' must be a single value")
  if(any(isUnknown(x, unknown=unknown)) && !force)
    stop(sprintf("'x' already has value %s", dQuote(unknown)))
  classX <- class(x)[1]
  classUnk <- class(unknown)[1]
  if(classX != classUnk) {
    tmp <- c("integer", "numeric")
    if(!(classX %in% tmp && classUnk %in% tmp)) {
      warning(sprintf("'unknown' should be %s for %s 'x' - will try to coerce",
                      dQuote(classX), dQuote(classX)), call.=call.)
    }
    unknown <- do.call(paste("as.", classX, sep=""), args=list(unknown))
  }
  x[is.na(x)] <- unknown
  x
}

NAToUnknown.factor <- function(x, unknown, force=FALSE, call.=FALSE, ...)
{
  if(length(unknown) != 1)
    stop("'unknown' must be a single value")
  if(any(isUnknown(x, unknown=unknown))) {
    if(!force) stop(sprintf("'x' already has level %s", dQuote(unknown)))
  } else {
    mapLevels(x) <- c(mapLevels(x, codes=FALSE),
                      mapLevels(as.character(unknown), codes=FALSE))
  }
  x[is.na(x)] <- unknown
  if(!force)
    warning(sprintf("new level is introduced: %s", unknown), call.=call.)
  x
}

NAToUnknown.list <- function(x, unknown, force=FALSE, call.=FALSE, ...)
{
  unknown <- .unknownList(x=x, unknown=unknown)
  x <- mapply(FUN="NAToUnknown", x=x, unknown=unknown, force=force,
              call.=call., SIMPLIFY=FALSE)
  x
}

NAToUnknown.data.frame <- function(x, unknown, force=FALSE, call.=FALSE, ...)
{
  x[] <- NAToUnknown.list(x=x, unknown=unknown, force=force, call.=call.)
  x
}

### }}}
### {{{ .unknownList
###------------------------------------------------------------------------

.unknownList <- function(x, unknown)
{
  ## --- Setup ---

  n <- length(x)
  unkN <- length(unknown)
  namesX <- names(x)
  namesXNullTest <- is.null(namesX)
  unkNames <- names(unknown)
  unkNamesNullTest <- is.null(unkNames)
  defInNames <- ".default" %in% unkNames
  defInd <- unkNames %in% ".default"
  def <- unknown[defInd]

  if(defInNames) { ## Remove default
    unkN <- unkN - 1
    unkNames <- unkNames[!defInd]
    unknown <- unknown[!defInd]
  }

  if(!namesXNullTest) { ## Check for nonexistent name
    test <- !(unkNames %in% namesX)
    if(any(test)) stop(sprintf("name(s) %s not in names of 'x'",
                       paste(sQuote(unkNames[test]), collapse=" ")))
  }

  ## --- Recycle ---

  if(unkN < n) {
    if(unkNamesNullTest | defInNames) {
      if(defInNames) { # handling .default
        names(def) <- NULL
        unknownDef <- rep(def, length=(n - unkN))
        names(unknownDef) <- namesX[!(namesX %in% unkNames)]
        unknown <- c(unknownDef, unknown)
      } else {
        unknownDef <- unknown
        unknown <- rep(unknownDef, length=n)
      }
    } else {
      stop("can not propely recycle named 'unknown'")
    }
  }

  ## --- Names ---

  if(!namesXNullTest) { ## no need if namesX NULL
    if(unkNamesNullTest) { ## missing unkNames
      names(unknown) <- namesX
    } else {                ## unkNames known
      unknown <- unknown[match(namesX, names(unknown))]
    }
  }

  unknown
}

### }}}
### {{{ Dear Emacs
### Local variables:
### folded-file: t
### End:
### }}}

###------------------------------------------------------------------------
### unknown.R ends here
