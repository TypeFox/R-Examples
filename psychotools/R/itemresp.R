## itemresp() is the class constructor:
##   - data: matrix/data.frame of n subjects (rows) and k items (columns)
##   - mscale: response categories / measurement scale for items (integer or character)
##   - labels: item labels
##   - names: subject names, IDs

itemresp <- function(data, mscale = NULL, labels = NULL, names = NULL)
{
  ## currently only used internally
  ncat <- NULL

  ## transform input to matrix/data.frame if provided as a vector
  if(is.vector(data)) data <- matrix(data, ncol = 1, dimnames = list(names(data), "item"))
  if(is.factor(data)) data <- data.frame(item = data, row.names = names(data))
  k <- ncol(data)

  ## check if numeric elements can be coerced to integer without loss
  num <- if(is.data.frame(data)) which(!sapply(data, is.factor)) else 1L:ncol(data)
  for (i in num) {
    if(!isTRUE(all.equal(as.integer(data[,i]), data[,i]))) stop(sprintf("item %i cannot be coerced to integer", num[i]))
  }

  ## item/subject labels
  if(!is.null(names)) rownames(data) <- names
  if(!is.null(labels)) colnames(data) <- labels
  if(is.null(colnames(data))) colnames(data) <- paste("item",
    formatC(1L:k, flag = "0", width = floor(log10(k)) + 1), sep = "")
  nam <- colnames(data)

  ## expand ncat to named list of length k
  ## (collapse to vector at the end)  
  ncat0 <- structure(vector(mode = "list", length = k), names = nam)
  if(is.null(ncat)) {
    ncat <- ncat0
  } else {
    ## if vector, turn into list of same length
    if(!is.list(ncat)) ncat <- as.list(ncat)

    ## if named, use names
    if(is.null(names(ncat))) {
      if(length(ncat) == 1L) ncat <- rep.int(ncat, k)
      if(length(ncat) != k) stop("length of ncat does not match number of items")
      ncat <- structure(ncat, names = nam)
    } else {
      if(!all(names(ncat) %in% c(nam, ""))) stop("names of ncat do not match names of items")
      if(any(table(names(ncat)) > 1L)) stop("names of ncat are not unique")
      if(any(names(ncat) == "")) {
        wi <- which(names(ncat) == "")
        ncat0[] <- ncat[wi]
	ncat <- ncat[-wi]
      }
      ncat0[names(ncat)] <- ncat
      ncat <- ncat0
    }
  }

  ## expand mscale to named list of length k
  mscale0 <- structure(vector(mode = "list", length = k), names = nam)
  if(is.null(mscale)) {
    mscale <- mscale0
  } else {
    ## if vector, turn into list of length 1
    if(!is.list(mscale)) mscale <- list(mscale)

    ## if named, use names
    if(is.null(names(mscale))) {
      if(length(mscale) == 1L) mscale <- rep.int(mscale, k)
      if(length(mscale) != k) stop("length of mscale does not match number of items")
      mscale <- structure(mscale, names = nam)
    } else {
      if(!all(names(mscale) %in% c(nam, ""))) stop("names of mscale do not match names of items")
      if(any(table(names(mscale)) > 1L)) stop("names of mscale are not unique")
      if(any(names(mscale) == "")) {
        wi <- which(names(mscale) == "")
        mscale0[] <- mscale[wi]
	mscale <- mscale[-wi]
      }
      mscale0[names(mscale)] <- mscale
      mscale <- mscale0
    }    
  }

  ## work though individual items
  for (i in 1L:k) {
    ## assure that ncat is either NULL or integer
    if(!is.null(ncat[[i]])) ncat[[i]] <- as.integer(ncat[[i]])

    ## bring ncat and mscale into sync
    ## (either both NULL or both specified)
    if(!is.null(ncat[[i]])) {
      if(is.null(mscale[[i]])) {
        mscale[[i]] <- 0L:(ncat[[i]] - 1L)
      } else {
        if(ncat[[i]] != length(mscale[[i]]))
	  stop(sprintf("length of mscale and number of categories does not match for item %i", i))
      }
    } else if(!is.null(mscale[[i]])) {
      ncat[[i]] <- length(mscale[[i]])
    }

    ## handle column (factor vs. numeric)
    if(is.factor(data[,i])) {
      if(is.null(mscale[[i]])) {
        mscale[[i]] <- levels(data[,i])
	ncat[[i]] <- length(mscale[[i]])
      } else {
        if(ncat[[i]] != length(levels(data[,i]))) stop(
	  sprintf("length of mscale and factor levels for item %i do not match", i))
	if(any(mscale[[i]] != levels(data[,i]))) warning(
	  sprintf("labels and factor levels for item %i do no match: mscale are used", i))
      }
      data[,i] <- as.numeric(data[,i]) - 1L
    } else {
      if(is.null(ncat[[i]])) {
        if(min(data[,i], na.rm = TRUE) < 0) {
          mscale[[i]] <- seq(-max(abs(data[,i]), na.rm = TRUE), max(abs(data[,i]), na.rm = TRUE), by = 1L)
        } else {
          mscale[[i]] <- seq(min(data[,i], na.rm = TRUE), max(data[,i], na.rm = TRUE), by = 1L)
        }
        ncat[[i]] <- length(mscale[[i]])      
        data[,i] <- data[,i] - min(mscale[[i]], na.rm = TRUE)
      } else {
        data[,i] <- as.numeric(factor(data[,i], levels = mscale[[i]])) - 1L
      }
    }

    ## check that ncat is greater 1
    if(ncat[[i]] <= 1L) stop(sprintf("item %i does not have more than one category", i))
  }
  
  ## coerce data to integer matrix (keep dim and dimnames)
  data <- structure(as.integer(as.matrix(data)), .Dim = dim(data), .Dimnames = dimnames(data))

  ## simplify mscale classes to integer if possible
  lclass <- sapply(mscale, class)
  if(any(lclass != "integer")) for(i in which(lclass != "integer")) {
    labi <- suppressWarnings(as.integer(mscale[[i]]))
    ok <- if(lclass[i] == "character") !any(is.na(labi)) else all(labi == mscale[[i]])
    if(ok) mscale[[i]] <- labi
  }
  lclass <- sapply(mscale, class)
  if(any(lclass != "integer")) mscale <- lapply(mscale, as.character)  
  
  ## attach attributes, class and return
  ## attr(data, "ncat") <- ncat
  attr(data, "mscale") <- mscale
  class(data) <- "itemresp"
  return(data)
}


##########

## type checking
is.itemresp <- function (x) inherits(x, "itemresp")
  
## subsetting/indexing
"[.itemresp" <- function(x, i, j, ...) {
  ms <- attr(x, "mscale")
  x <- unclass(x)[i, j, drop = FALSE]
  attr(x, "mscale") <- if(missing(j)) ms else ms[j]
  structure(x, class = "itemresp")
} 

subset.itemresp <- function(x, items = NULL, subjects = NULL, ...) {
  if(is.null(items)) items <- 1L:ncol(x)
  if(is.null(subjects)) subjects <- 1L:nrow(x)
  x[i = subjects, j = items]
}

length.itemresp <- function(x) nrow(unclass(x))

c.itemresp <- function(...)
{
  args <- list(...)

  check_list <- function(x) {
    if(length(x) < 2L) return(TRUE)
    x1 <- x[[1L]]
    all(sapply(2L:length(x), function(i) identical(x[[i]], x1)))
  }
  if(!check_list(lapply(args, colnames))) stop("objects have different item labels")
  if(!check_list(lapply(args, attr, "mscale"))) stop("objects have different mscales")

  ms <- attr(args[[1L]], "mscale")
  rval <- do.call("rbind", lapply(args, function(x) as.matrix(x)))
  attr(rval, "mscale") <- ms
  structure(rval, class = "itemresp")
}

merge.itemresp <- function(x, y, all = FALSE, ...)
{
  ## check item labels
  if(any(colnames(x) %in% colnames(y))) {
    colnames(x) <- paste(deparse(substitute(x)), colnames(x), sep = ".")
    colnames(y) <- paste(deparse(substitute(y)), colnames(y), sep = ".")
  }

  ## check subjects
  xnam <- rownames(x)
  ynam <- rownames(y)
  nonames <- FALSE
  if(is.null(xnam) | is.null(ynam)) {
    stopifnot(nrow(x) == nrow(y))
    if(is.null(xnam)) xnam <- ynam
    if(is.null(xnam)) {
      xnam <- as.character(1:nrow(x))
      nonames <- TRUE
    }
    if(is.null(ynam)) ynam <- xnam
  }
  
  ## coerce to data frame and merge along auxiliary name column
  x <- as.list(x, df = TRUE)
  y <- as.list(y, df = TRUE)
  x[["(.names)"]] <- xnam
  y[["(.names)"]] <- ynam
  xy <- merge(x, y, by = "(.names)", all = all, ...)
  rownames(xy) <- xy[["(.names)"]]
  xy[["(.names)"]] <- NULL
  xy <- xy[intersect(c(xnam, ynam), rownames(xy)), , drop = FALSE]

  ## transform to itemresp and return
  rval <- itemresp(xy)
  if(nonames) names(rval) <- NULL
  return(rval)
}

rep.itemresp <- function(x, ...) {
  ix <- if(length(x) > 0L) 1L:length(x) else 0L
  x[rep(ix, ...)]
}

xtfrm.itemresp <- function(x) {
  if(length(x) > 0L) 1L:length(x) else 0L
}


##########
  

## extract/replace mscale
mscale.itemresp <- function(object, items = NULL, simplify = NULL, ...)
{
  ## extract full list
  ms <- attr(object, "mscale")

  ## items for subsetting
  if(is.null(items)) items <- colnames(object)
  if(!is.character(items)) items <- colnames(object)[items]

  ## vector if single item and !simplify
  if(length(items) == 1L) {
    if(identical(simplify, FALSE)) return(ms[items]) else return(ms[[items]])
  }
  
  ## list otherwise
  rval <- ms[items]
  ## or try to simplify
  if(!identical(simplify, FALSE) && all(sapply(seq_along(rval), function(i) identical(rval[[1L]], rval[[i]])))) {
    rval <- as.vector(rval[[1L]])
  } else {
    if(identical(simplify, TRUE)) warning("mscale could not be simplified, list returned")
  }

  return(rval)
}
levels.itemresp <- function(x, ...) mscale(x, ...)

"mscale<-.itemresp" <- function(object, value) {
  k <- ncol(object)
  nam <- colnames(object)
  
  ## handle "new" mscale value
  value0 <- attr(object, "mscale")
  if(!is.list(value)) value <- list(value)
  if(is.null(names(value))) {
    if(length(value) == 1L) value <- rep.int(value, k)
    if(length(value) != k) stop("length of mscale does not match number of items")
    value <- structure(value, names = nam)
  } else {
    if(!all(names(value) %in% c(nam, ""))) stop("names of mscale do not match names of items")
    if(any(table(names(value)) > 1L)) stop("names of mscale are not unique")
    if(any(names(value) == "")) {
      wi <- which(names(value) == "")
      value0[] <- value[wi]
      value <- value[-wi]
    }
    value0[names(value)] <- value
    value <- value0
  }

  ## replace mscale (if correct length)
  ms <- attr(object, "mscale")
  ok <- sapply(ms, length) == sapply(value, length)
  if(!all(ok)) warning(sprintf("mscale length incorrect for item %s; not replaced", paste(which(!ok), collapse = ", ")))
  ms[ok] <- value[ok]

  ## clean up integer vs. character, simplify to integer if possible
  lclass <- sapply(ms, class)
  if(any(lclass != "integer")) for(i in which(lclass != "integer")) {
    labi <- suppressWarnings(as.integer(ms[[i]]))
    ok <- if(lclass[i] == "character") !any(is.na(labi)) else all(labi == ms[[i]])
    if(ok) ms[[i]] <- labi
  }
  lclass <- sapply(ms, class)
  if(any(lclass != "integer")) ms <- lapply(ms, as.character)  
  
  ## check for duplicated values (for collapsing responses)
  dup <- sapply(ms, function(x) any(duplicated(x)))
  if(any(dup)) {
    object <- unclass(object)  
    for(i in which(dup)) {
      object[, i] <- as.integer(factor(ms[[i]], levels = unique(ms[[i]])))[object[, i] + 1L] - 1L
      ms[[i]] <- unique(ms[[i]])
    }
    class(object) <- "itemresp"
  }
  
  ## actually replace attribute and return
  attr(object, "mscale") <- ms
  return(object)
}

## names = subject names = stored as row names internally
names.itemresp <- function(x) rownames(x)

"names<-.itemresp" <- function(x, value) {
  rownames(x) <- value
  return(x)
}

## labels = item names = stored as column names internally
labels.itemresp <- function(object, ...) colnames(unclass(object))

"labels<-.itemresp" <- function(object, value) {
  colnames(object) <- names(attr(object, "mscale")) <- value
  return(object)
}


##########
  

## Methods: format() handles character formatting
##   utilized in print() and as.character() method.
format.itemresp <- function(x, sep = c(",", ":"), brackets = TRUE,
  abbreviate = NULL, mscale = TRUE, labels = FALSE, width = getOption("width") - 7L, ...)
{
  ## process brackets
  if(is.null(brackets)) brackets <- TRUE
  if(is.logical(brackets)) brackets <- if(brackets) c("{", "}") else ""
  brackets <- rep(as.character(brackets), length.out = 2L)

  ## process separators
  sep <- rep(as.character(sep), length.out = 2L)

  ## internal building blocks
  ms <- attr(x, "mscale")
  x <- as.matrix(x)

  ## abbreviate
  ms <- lapply(ms, as.character)
  maxchar <- max(unlist(lapply(ms, nchar)))
  if(is.null(abbreviate)) abbreviate <- mscale & maxchar > 3L
  if(!identical(abbreviate, FALSE)) {
    minlength <- if(is.logical(abbreviate)) {
      as.numeric(cut(maxchar, c(-Inf, 1.5, 4.5, 7.5, Inf)))
    } else {
      abbreviate
    }
    ms <- lapply(ms, abbreviate, minlength = minlength)
  }
  
  ## coerce to character
  rval <- lapply(seq_along(ms), function(i) paste(
    if(labels) paste(colnames(x)[i], sep[2L], sep = "") else "",
    if(mscale) ms[[i]][x[, i] + 1L] else as.character(x[, i]),
    sep = ""
  ))
  rval <- do.call("cbind", rval)
  
  ## collapse over items
  rval <- apply(rval, 1L, paste, collapse = sep[1L])
  
  ## check width
  if(is.null(width)) width <- TRUE
  if(is.logical(width)) width <- if(width) getOption("width") - 7L else Inf
  width <- width - sum(nchar(brackets))
  wi <- nchar(rval) > width
  if(any(wi)) {
    rval[wi] <- paste(substr(rval[wi], 1L, width - 3L), "...", sep = "")
  }
  
  ## add brackets
  rval <- paste(brackets[1L], rval, brackets[2L], sep = "")

  ## assure names (if any)
  names(rval) <- rownames(x)

  return(rval)
}

print.itemresp <- function(x, quote = FALSE, ...)
{
  print(format(x, ...), quote = quote)
  invisible(x)
}


##########
  

## coercion functions
as.character.itemresp <- function(x, ...) {
  format(x, mscale = TRUE, labels = FALSE, width = FALSE, ...)
}

as.integer.itemresp <- 
as.matrix.itemresp <- function(x, ...) {
  rval <- unclass(x)
  attr(rval, "mscale") <- NULL
  class(rval) <- "matrix"
  return(rval)  
}

as.double.itemresp <- function(x, ...) {
  rval <- as.matrix(x, ...)
  rval[] <- as.double(rval)
  return(rval)
}

as.data.frame.itemresp <- function(x, row.names = NULL, optional = FALSE, ...,
  nm = paste(deparse(substitute(x), width.cutoff = 500L), collapse = " "))
{
  force(nm)
  nrows <- length(x)
  if(is.null(row.names)) {
    if(nrows == 0L) { 
      row.names <- character()
    } else {
      if(length(row.names <- names(x)) == nrows && !anyDuplicated(row.names)) {
      } else {
        row.names <- .set_row_names(nrows)
      }
    }
  }
  if(!is.null(names(x))) names(x) <- NULL
  value <- list(x)
  if(!optional) names(value) <- nm
  attr(value, "row.names") <- row.names
  class(value) <- "data.frame"
  value
}

as.list.itemresp <- function(x, items = NULL, mscale = TRUE, df = FALSE, ...)
{
  ## building blocks
  ms <- attr(x, "mscale")
  x <- unclass(x)

  ## select desired items
  if(!is.null(items)) {
    ms <- ms[items]
    x <- x[, items, drop = FALSE]
  }

  ## turn matrix into list of factors
  rval <- lapply(1L:ncol(x), function(i) factor(
    x[, i],
    levels = 0L:(length(ms[[i]]) - 1L),
    labels = if(mscale) ms[[i]] else 0L:(length(ms[[i]]) - 1L)
  ))
  names(rval) <- colnames(x)
  
  ## coerce to data.frame if requested
  if(df) rval <- as.data.frame(rval)
  
  ## return
  return(rval)
}


##########
  

## just shows structure (no data),
str.itemresp <- function(object, width = getOption("width") - 7L, ...)
{
  rval <- c(
    sprintf(" Item response data from %d subjects for %d items.\n",
      nrow(unclass(object)),
      ncol(unclass(object))),
    paste(labels(object), ": ", sapply(attr(object, "mscale"), paste, collapse = ", "), "\n", sep = "")
  )
  
  if(is.null(width)) width <- TRUE
  if(is.logical(width)) width <- if(width) getOption("width") - 7L else Inf
  ok <- nchar(rval) <= width
  if(any(!ok)) rval[!ok] <- paste(substr(rval[!ok], 1, width - 3L), "...\n", sep = "")
  
  cat(rval)
  invisible(NULL)
}

## summary() computes frequency table
summary.itemresp <- function(object, items = NULL, abbreviate = FALSE, mscale = TRUE, simplify = TRUE, sep = " ", ...)
{
  ## frequency table for each item
  rval <- as.list(object, items = items, mscale = mscale)
  rval <- lapply(rval, table)  
  if(!simplify) return(rval)
  
  ## abbreviate
  ms <- attr(object, "mscale")
  ms <- lapply(ms, as.character)
  maxchar <- max(unlist(lapply(ms, nchar)))
  if(is.null(abbreviate)) abbreviate <- mscale & maxchar > 3L
  if(!identical(abbreviate, FALSE)) {
    minlength <- if(is.logical(abbreviate)) {
      as.numeric(cut(max(unlist(lapply(ms, nchar))), c(-Inf, 1.5, 4.5, 7.5, Inf)))
    } else {
      abbreviate
    }
    ms <- lapply(ms, abbreviate, minlength = minlength)
  }

  ## collapse to table
  if(!is.null(items)) ms <- ms[items]
  ncat <- sapply(ms, length)
  if(!all(ncat == ncat[1L])) for(i in which(ncat < max(ncat))) {
      rval[[i]] <- c(rval[[i]], rep.int(NA_integer_, max(ncat) - ncat[[i]]))
  }
  rval <- do.call("rbind", rval)
  if(mscale & all(sapply(seq_along(ms), function(i) identical(ms[[1L]], ms[[i]])))) {
    colnames(rval) <- ms[[1L]]
  } else {
    colnames(rval) <- 0:(max(ncat) - 1L)
    if(mscale) rownames(rval) <- paste(rownames(rval),
      paste("(", sapply(ms, paste, collapse = ","), ")", sep = ""),
      sep = sep)
  }
  return(rval)
}

plot.itemresp <- function(x, xlab = "", ylab = "", items = NULL, abbreviate = FALSE, mscale = TRUE, sep = "\n",
  off = 2, axes = TRUE, names = TRUE, srt = 45, adj = c(1.1, 1.1), ...)
{
  ## summarize item response data
  tab <- summary(x, items = items, abbreviate = abbreviate, mscale = mscale, simplify = TRUE, sep = sep)
  tab[is.na(tab)] <- 0
  n <- NROW(tab)

  ## process labeling option for item names
  if (isTRUE(names)) nms <- rownames(tab)
  if (is.character(names)) {
    nms <- names
    names <- TRUE
  }
  if(!names) {
    lab <- rep(NA, n)
    lab[c(1, n)] <- c(1, n)
    pr <- pretty(1:n)
    pr <- pr[pr > 1 & pr < n]
    lab[pr] <- pr    
    nms <- lab
  }
  lab <- rownames(tab)
  rownames(tab) <- rep.int("", n)

  ## call spineplot
  spineplot(tab, xlab = xlab, ylab = ylab, off = off, axes = axes, ...)

  ## x-axis labels
  ix <- rowSums(tab)/sum(tab)
  ix <- cumsum(ix + off/100) - off/100 - ix/2
  if(names) {
    text(ix, par("usr")[3], labels = nms, srt = srt, adj = adj, xpd = TRUE, cex = 0.9)
  } else {
    axis(1, at = ix, labels = nms)
  }

  ## return summary table invisibly
  rownames(tab) <- lab
  invisible(tab)
}


##########
  

## default NA handling
is.na.itemresp <- function(x) apply(is.na(unclass(x)), 1L, all)

