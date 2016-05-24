## util.R (part of the TRAMPR package)

## These are small utility functions.  I will document these, but
## probably won't expose them.

## Assert that a data.frame must contain columns
must.contain.cols <- function(d, cols) {
  i <- cols %in% names(d)
  if ( !all(i) )
    stop(sprintf("Required columns missing in %s: %s",
                 dQuote(deparse(substitute(d))),
                 paste(dQuote(cols[!i]), collapse=", ")))
  invisible(TRUE)
}

## This is like match() for data.frames.
## Given a data.frame 'x' and a table 'table', for every row in 'x',
## return the the row number where all columns in 'x' that are in
## 'table' match between 'x' and 'table'.
## x: a b c  table: a b   -> result:
##    1 2 9         1 2   ->         1
##    1 2 8         2 2   ->         1
##    2 2 7               ->         2
## This uses the same method as duplicated() for combining columns
## (which differs from the method used by tapply).  See notes in
## ?duplicated for pitfalls.
classify <- function(x, table, ...) {
  x.str <- do.call(paste, c(x[names(table)], sep="\r"))
  x.str[rowSums(is.na(x[names(table)])) > 0] <- NA
  table.str <- do.call(paste, c(table, sep="\r"))
  table.str[rowSums(is.na(table)) > 0] <- NA
  match(x.str, table.str, ...)
}

## Return the value of the element with the minimum _absolute_
## value, but including the original sign (-.1 is less than -1 or 1)
## If all values are NA, then return NA (rather than numeric(0), which
## which.min() defaults to).
absolute.min <- function(x)
  if ( all(is.na(x)) && length(x) ) NA else x[which.min(abs(x))]

## Read in a string, with a default, perhaps checking.
read.string <- function(prompt, default=NULL, check=TRUE) {
  repeat {
    s <- strstrip(readline(prompt))
    if ( nchar(s) == 0 ) {
      if ( missing(default) )
        next
      s <- default
      break
    }
    if ( check )
      break
    
    action <-
      strstrip(readline(sprintf("Recieved `%s'; accept [Y/n] ", s)))
    if ( tolower(substr(action, 1, 1)) != "n" )
      break
  }
  s
}

## Strip leading and trailing whitespace from a string
strstrip <- function(s)
  sub("^\\s*(?U)(.*)\\s*$", "\\1", s, perl=TRUE)

## Join two data.frames that share some (but perhaps not all) columns:
## row.names handling may change with R 2.4.0
rbind2 <- function(x, y) {
  x[setdiff(names(y), names(x))] <- NA
  y[setdiff(names(x), names(y))] <- NA
  rbind(x, y)
}

defactor <- function(x, warn=FALSE) {
  stopifnot(is.data.frame(x))
  i <- sapply(x, is.factor)
  if ( any(i) ) {
    if ( warn )
      warning("Converting ", sum(i), " columns from factors: ",
              paste(names(x)[i], collapse=", "))
    x[i] <- lapply(x[i], as.character)
  }
  x
}

## New defaults for read.csv:
read.csv.safe <- function(file, ...)
  read.csv(file, as.is=TRUE, strip.white=TRUE, na.strings=c("NA", ""),
           ...)

## Wrap labels so they fit within a space.  Probably a little more
## useful than just for TRAMPR.
labwrap <- function(x, width, indent=0, exdent=0, cex=1,
                    units="user") {
  if (!is.character(x)) 
    x <- as.character(x)
  indentString <- paste(rep.int(" ", indent), collapse = "")
  exdentString <- paste(rep.int(" ", exdent), collapse = "")
  spaceWidth <- strwidth(" ", units, cex)
  indentWidth <- strwidth(indentString, units, cex)
  exdentWidth <- strwidth(exdentString, units, cex)

  lenFirst <- width - indentWidth
  lenRest <- width - exdentWidth

  z <- strsplit(x, "[ \t\n]")
  y <- list()

  for ( i in seq_along(z) ) {
    words <- z[[i]]
    nwords <- length(words)
    wid <- strwidth(words, units, cex)

    idx <- rep.int(0, nwords)
    idx.word <- 1

    for ( line in seq_len(nwords) ) {
      lenAllowed <- if ( line == 1 ) lenFirst else lenRest
      lens <- cumsum(wid[idx.word:nwords] + spaceWidth)
      over <- lens > lenAllowed
      done <- !any(over)

      at <- if ( done ) length(lens) else max(2, which(over)[1]) - 1
      
      thisLine <- seq(idx.word, length=at)
      idx[thisLine] <- line

      if ( done )
        break
      else
        idx.word <- idx.word + at
    }

    ## Now construct the strings:
    s <- lapply(split(words, idx), paste, collapse=" ")
    indent.exdent <- rep(c(indentString, exdentString),
                         c(1, line-1))
    y[[i]] <- paste(indent.exdent, s, sep="")
  }
  y
}
