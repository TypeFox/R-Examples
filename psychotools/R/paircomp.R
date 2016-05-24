## paircomp() is the class constructor:
##   - data: integer coded matrix of n subjects (rows) and
##            choose(k, 2) paired comparisons (columns)
##   - labels: labels of k objects
##   - mscale: measurement scale for comparisons (centered around 0)
##   - ordered: logical. Is a:b different from b:a?
##   - covariates: data frame with k rows for object covariates

paircomp <- function(data, labels = NULL, mscale = NULL, ordered = FALSE, covariates = NULL)
{
  ## data should be matrix(-like)
  npc <- NCOL(data)
  if(npc == 1L | is.data.frame(data)) data <- as.matrix(data)

  ## compute unique non-NA values
  data_unique <- as.vector(na.omit(unique(as.vector(data))))

  ## compute number of subjects and objects
  nsubj <- NROW(data)
  nobj <- 0.5 + sqrt(0.25 + if(ordered) npc else 2 * npc)
  
  ## covariates should be data.frame(-like)
  if(!is.null(covariates)) covariates <- as.data.frame(covariates)
  
  ## sanity checks
  stopifnot(is.matrix(data), nobj >= 2, isTRUE(all.equal(nobj, round(nobj))),
    isTRUE(all.equal(data_unique, round(data_unique))))
  if(!is.null(covariates)) stopifnot(nrow(covariates) == nobj)

  ## coerce to integer and set dimnames
  cnam <- which(upper.tri(diag(nobj)), arr.ind = TRUE)
  if(ordered) cnam <- rbind(cnam, cnam[,2:1])
  cnam <- apply(cnam, 1, paste, collapse = ":")
  data <- structure(as.integer(data), .Dim = c(nsubj, npc), .Dimnames = list(rownames(data), cnam))
  data_unique <- as.integer(data_unique)

  ## process labels
  if(is.null(labels)) labels <- make.unique(rep(letters, length.out = nobj), sep = "") else {
    if(length(labels) != nobj) stop("length of labels does not match number of objects")
  }

  ## process mscale
  if(is.null(mscale)) {
    mscale <- if(length(data_unique) < 1) c(-1, 1) else {
      mscale <- max(abs(data_unique))
      if(!(mscale > 0)) mscale <- 1
      mscale <- (-mscale):mscale
      if(!(0 %in% data_unique)) mscale <- mscale[-(length(mscale)+1)/2]
      mscale
      }
  } else {
    if(!all(data_unique %in% mscale)) stop("mscale does not match data")
    mscale <- as.integer(sort(mscale))
    if(max(abs(mscale)) <= 0) stop("mscale needs to have non-zero elements")
    if(abs(head(mscale, 1)) != tail(mscale, 1)) stop("mscale must by symmetric")
  }

  ## process covariates
  if(!is.null(covariates)) rownames(covariates) <- seq_along(labels)

  rval <- data
  class(rval) <- "paircomp"
  attributes(rval) <- c(attributes(rval), list(
    labels = labels,
    mscale = mscale,
    ordered = ordered,
    covariates = covariates))
  class(rval) <- "paircomp"
  return(rval)
}



## Methods: format() handles character formatting
##   utilized in print() and as.character() method.

format.paircomp <- function(x, sep = ", ", brackets = TRUE,
  abbreviate = NULL, width = getOption("width") - 7, ...)
{
  ## process brackets
  if(is.null(brackets)) brackets <- TRUE
  if(is.logical(brackets)) brackets <- if(brackets) c("{", "}") else ""
  brackets <- rep(as.character(brackets), length.out = 2)

  ## set up comparison symbols
  mscale <- as.integer(sort(unique(c(mscale(x), 0))))
  mscale_max <- tail(mscale, 1)
  mscale_symbol <- if(mscale_max > 2L) {
    paste(ifelse(abs(mscale) > 1L, abs(mscale), ""), c("<", "=", ">")[sign(mscale) + 2L], sep = "")
  } else if(mscale_max == 2L) {
    c("<<", "<", "=", ">", ">>")[mscale + 3L]
  } else {
    c("<", "=", ">")[mscale + 2L]
  }

  ## process labels
  lab <- labels(x)
  ## abbreviate
  if(is.null(abbreviate)) abbreviate <- TRUE
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(lab))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  lab <- abbreviate(lab, abbreviate)  
  ## expand
  ix <- which(upper.tri(diag(length(lab))), arr.ind = TRUE)
  if(!attr(x, "ordered")) {
    lab1 <- lab[ix[,1]]
    lab2 <- lab[ix[,2]]
  } else {
    lab1 <- c(lab[ix[,1]], lab[ix[,2]])
    lab2 <- c(lab[ix[,2]], lab[ix[,1]])
  }
  
  pc <- as.matrix(x)
  pc <- pc + mscale_max + 1
  pc <- apply(pc, 1, function(x) paste(
    ifelse(is.na(x), "NA", paste(lab1, mscale_symbol[x], lab2, sep = " ")),
    collapse = sep))
  ## Handle NAs differently? Maybe via "a ? b"?

  ## check width
  if(is.null(width)) width <- TRUE
  if(is.logical(width)) width <- if(width) getOption("width") - 7 else Inf
  width <- width - sum(nchar(brackets))
  wi <- nchar(pc) > width
  if(any(wi)) {
    pc[wi] <- paste(substr(pc[wi], 1, width - 3), "...", sep = "")
  }
  
  ## add brackets
  pc <- paste(brackets[1], pc, brackets[2], sep = "")

  ## assure names (if any)
  names(pc) <- rownames(unclass(x))

  return(pc)
}

print.paircomp <- function(x, quote = FALSE, ...)
{
  print(format(x, ...), quote = quote)
  invisible(x)
}



## Basic summaries:
## str() just shows structure (no data),
## summary() shows frequency table for each comparison.

str.paircomp <- function(object, width = getOption("width") - 7, ...)
{
  rval <- sprintf(" Paired comparisons from %d subjects for %d objects: %s.\n",
    nrow(unclass(object)), length(labels(object)), paste(labels(object), collapse = ", "))
  
  if(is.null(width)) width <- TRUE
  if(is.logical(width)) width <- if(width) getOption("width") - 7 else Inf
  if(nchar(rval) > width) rval <- paste(substr(rval, 1, width - 3), "...\n", sep = "")
  
  cat(rval)
  invisible(NULL)
}

summary.paircomp <- function(object, abbreviate = FALSE, decreasing = TRUE, pcmatrix = FALSE, weights = NULL, ...)
{
  ## data
  dat <- as.matrix(object)
  
  ## process labels
  lab <- labels(object)
  ## abbreviate
  if(is.null(abbreviate)) abbreviate <- TRUE
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(lab))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  lab <- abbreviate(lab, abbreviate)  

  ## rownames
  ix <- which(upper.tri(diag(length(lab))), arr.ind = TRUE)
  if(!attr(object, "ordered")) {
    lab1 <- lab[ix[,1]]
    lab2 <- lab[ix[,2]]
  } else {
    lab1 <- c(lab[ix[,1]], lab[ix[,2]])
    lab2 <- c(lab[ix[,2]], lab[ix[,1]])
  }
  rnam <- paste(format(lab1), ":", format(lab2))

  ## colnames
  mscale <- mscale(object)
  mscale_max <- tail(mscale, 1)
  cnam <- if(mscale_max > 2L) {
    paste(ifelse(abs(mscale) > 1L, abs(mscale), ""), c("<", "=", ">")[sign(mscale) + 2L], sep = "")
  } else if(mscale_max == 2L) {
    c("<<", "<", "=", ">", ">>")[mscale + 3L]
  } else {
    c("<", "=", ">")[mscale + 2L]
  }
  if(decreasing) cnam <- rev(cnam)
  if(any(is.na(dat))) cnam <- c(cnam, "NA's")

  rval <- if(is.null(weights)) {
    t(apply(dat, 2, function(x) table(factor(x, levels = mscale(object)))))
  } else {
    weights <- rep(weights, length.out = nrow(dat))
    t(apply(dat, 2, function(x) xtabs(weights ~ factor(x, levels = mscale(object)))))
  }
  if(decreasing) rval <- rval[, ncol(rval):1, drop = FALSE]
  if(any(is.na(dat))) rval <- cbind(rval, apply(dat, 2, function(x) sum(is.na(x))))
  dimnames(rval) <- list(rnam, cnam)

  ## return paired-comparison matrix (only for unordered binary paircomp's)
  ## FIXME: better name for argument / standalone method?
  ## FIXME: table() and xtabs() are both not generic...
  if(pcmatrix & length(mscale) == 2 & !attr(object, "ordered")) {
    mat <- matrix(0, ncol = length(lab), nrow = length(lab))
    rownames(mat) <- colnames(mat) <- lab
    names(dimnames(mat)) <- cnam[1:2]
    mat[upper.tri(mat)] <- rval[,1]
    mat <- t(mat)
    mat[upper.tri(mat)] <- rval[,2]
    rval <- t(mat)
  }
  if(pcmatrix & length(mscale) == 2 & attr(object, "ordered")) {
    ordarr <- array(0, c(length(lab), length(lab), 2))

    mat <- matrix(0, ncol = length(lab), nrow = length(lab))
    mat[upper.tri(mat)] <- rval[1:(nrow(rval)/2), 1]
    mat <- t(mat)
    mat[upper.tri(mat)] <- rval[1:(nrow(rval)/2), 2]
    ordarr[,,1] <- t(mat)
    mat[upper.tri(mat)] <- rval[(nrow(rval)/2 + 1):nrow(rval), 2]
    mat <- t(mat)
    mat[upper.tri(mat)] <- rval[(nrow(rval)/2 + 1):nrow(rval), 1]
    ordarr[,,2] <- t(mat)

    dimnames(ordarr) <- list(lab, lab, c("1", "2"))
    names(dimnames(ordarr)) <- c(cnam[1:2], "order")
    rval <- ordarr
  }
  
  rval
}



## extract and combine observations

length.paircomp <- function(x) nrow(unclass(x))

"[.paircomp" <- function(x, i, ...) {
  xattr <- attributes(x)[c("labels", "mscale", "ordered", "covariates")]
  x <- unclass(x)[i,,drop=FALSE]
  attributes(x) <- c(attributes(x), xattr)
  structure(x, class = "paircomp")
}

c.paircomp <- function(...)
{
  args <- list(...)

  check_list <- function(x) {
    if(length(x) < 2) return(TRUE)
    x1 <- x[[1]]
    all(sapply(2:length(x), function(i) identical(x[[i]], x1)))
  }
  if(!check_list(lapply(args, attr, "labels"))) stop("objects have different labels")
  if(!check_list(lapply(args, attr, "mscale"))) stop("objects have different mscales")
  if(!check_list(lapply(args, attr, "ordered"))) stop("objects are differently ordered")
  if(!check_list(lapply(args, attr, "covariates"))) stop("objects have different covariates")

  rval <- args[[1]]
  xattr <- attributes(rval)[c("labels", "mscale", "ordered", "covariates")]
  rval <- do.call("rbind", lapply(args, function(x) as.matrix(x)))
  attributes(rval) <- c(attributes(rval), xattr)
  structure(rval, class = "paircomp")
}

rep.paircomp <- function(x, ...) {
  ix <- if(length(x) > 0) 1:length(x) else 0
  x[rep(ix, ...)]
}

xtfrm.paircomp <- function(x) {
  if(length(x) > 0) 1:length(x) else 0
}

## mscale(): new generic with extractor method for paircomp

mscale <- function(object, ...) UseMethod("mscale")

mscale.paircomp <- function(object, ...) attr(object, "mscale")

"mscale<-" <- function(object, value) UseMethod("mscale<-")

"mscale<-.paircomp" <- function(object, value) {
  ms <- attr(object, "mscale")
  val <- as.integer(sort(value))
  if(!isTRUE(all.equal(val, value))) warning("mscale sorted to be in increasing order")
  stopifnot(all(val %in% ms) | all(ms %in% val) | length(ms) == length(value))
  if(max(abs(val)) <= 0) stop("mscale needs to have non-zero elements")
  if(abs(head(val, 1)) != tail(val, 1)) stop("mscale must by symmetric")

  if(length(ms) == length(value)) {
    object[] <- val[as.vector(object) - min(ms) + 1]
  } else if(all(val %in% ms)) {
    object[!(as.vector(object) %in% val)] <- NA
  }

  attr(object, "mscale") <- unique(val)
  return(object)
}


## covariates(): new generic with extractor method for paircomp

covariates <- function(object, ...) UseMethod("covariates")

"covariates<-" <- function(object, value) UseMethod("covariates<-")

covariates.paircomp <- function(object, ...) {
  if(is.null(attr(object, "covariates"))) return(NULL)
  rval <- attr(object, "covariates")
  rownames(rval) <- labels(object)
  return(rval)
}

"covariates<-.paircomp" <- function(object, value) {
  if(!is.data.frame(value)) {
    dval <- as.data.frame(value)
    if(ncol(dval) == 1) names(dval) <- paste(deparse(substitute(value), width.cutoff = 500), collapse = " ")
    value <- dval
  }
  if(nrow(value) != length(labels(object))) stop("Number of rows in covariates does not match number of objects.")
  rownames(value) <- 1:length(labels(object))
  attr(object, "covariates") <- value
  return(object)
}


## names() queries/sets subject names,
## labels() queries/sets object names, (new generic labels<-)
## levels() is an alias for labels() (but discouraged to avoid confusion with mscale)

names.paircomp <- function(x) rownames(unclass(x))

"names<-.paircomp" <- function(x, value) {
  x <- unclass(x)
  rownames(x) <- value
  structure(x, class = "paircomp")
}

"labels<-" <- function(object, value) UseMethod("labels<-")
  
labels.paircomp <- function(object, ...) attr(object, "labels")

"labels<-.paircomp" <- function(object, value) {
  if(!(length(value) == length(attr(object, "labels")))) stop("length of labels does not match")
  attr(object, "labels") <- value
  object
}

## The idea was to have levels() synonymous with labels(),
## but there are situations where R then concludes that this
## is in fact a factor, argh...
##
## levels.paircomp <- function(x, ...) attr(x, "labels")
## 
## "levels<-.paircomp" <- function(x, value) {
##   labels(x) <- value
##   x
## }



## reorder() method enables selection of subsets
## of objects to be compared and/or re-ordering of
## objects.
## subset() currently just calls reorder(), in addition
## selection of subsets of subjects would be nice to have.

reorder.paircomp <- function(x, labels, ...)
{
  xlab <- labels(x)
  if(missing(labels)) labels <- xlab
  if(is.character(labels)) labels <- sapply(labels, match.arg, choices = xlab)
  if(is.character(labels)) labels <- match(labels, xlab)
  if(length(labels) < 2) stop("Need to compare at least two objects.")
  
  ## select relevant columns
  nlab <- labels
  ix <- which(upper.tri(diag(length(xlab))), arr.ind = TRUE)
  wi <- apply(ix, 1, function(z) all(z %in% nlab))
  ndat <- if(attr(x, "ordered")) as.matrix(x)[, c(wi, wi), drop = FALSE]
    else as.matrix(x)[, wi, drop = FALSE]
  
  ## re-order comparisons (if necessary)
  if(!identical(nlab, sort(nlab))) {
    if(!attr(x, "ordered")) {    
      ## set up indexes
      ix <- ix[wi,, drop = FALSE]
      nix <- which(upper.tri(diag(length(nlab))), arr.ind = TRUE)
      nix <- matrix(nlab[nix], ncol = 2)
      ## match
      wi <- apply(nix, 1, function(x) x[1] > x[2])
      nix[wi,] <- nix[wi, 2:1]
      onam <- apply(ix, 1, paste, collapse = ":")
      nnam <- apply(nix, 1, paste, collapse = ":")
      ord <- match(nnam, onam)
      ## reorder
      ndat <- ndat[, ord, drop = FALSE]
      ndat[, wi] <- -1 * ndat[, wi]
    } else {
      ## set up indexes
      ix <- ix[wi,, drop = FALSE]
      nix <- which(upper.tri(diag(length(nlab))), arr.ind = TRUE)
      nix <- matrix(nlab[nix], ncol = 2)
      ## match
      onam <- apply(rbind(ix, ix[,2:1]), 1, paste, collapse = ":")
      nnam <- apply(rbind(nix, nix[,2:1]), 1, paste, collapse = ":")
      ord <- match(nnam, onam)
      ## reorder
      ndat <- ndat[, ord, drop = FALSE]
    }    
  }

  ## call constructor to setup proper new paircomp object
  paircomp(ndat,
    labels = attr(x, "labels")[nlab], mscale = mscale(x),
    ordered = attr(x, "ordered"), covariates = attr(x, "covariates")[nlab,,drop=FALSE])
}

subset.paircomp <- function(x, subset, select, ...)
{
  if(!missing(subset)) stop("'subset' argument not yet implemented.")
  ## FIXME: This should enable the user to select all observations with, say:
  ## "a > b" or "a : b > -1" etc.

  reorder(x, labels = select, ...)
}


## coercion functions, all fairly straightforward

as.character.paircomp <- function(x, ...) {
  format(x, abbreviate = FALSE, width = FALSE, ...)
}

as.integer.paircomp <-
as.matrix.paircomp <- function(x, ...) {
  rval <- unclass(x)
  attr(rval, "labels") <- attr(rval, "mscale") <- attr(rval, "ordered") <- attr(rval, "covariates") <- NULL

  ## colnames
  lab <- attr(x, "labels")
  ix <- which(upper.tri(diag(length(lab))), arr.ind = TRUE)
  if(!attr(x, "ordered")) {
    lab1 <- lab[ix[,1]]
    lab2 <- lab[ix[,2]]
  } else {
    lab1 <- c(lab[ix[,1]], lab[ix[,2]])
    lab2 <- c(lab[ix[,2]], lab[ix[,1]])
  }
  colnames(rval) <- paste(lab1, ":", lab2, sep = "")
  return(rval)  
}

as.double.paircomp <- function(x, ...) {
  rval <- as.matrix(x, ...)
  rval[] <- as.double(rval)
  return(rval)
}

as.data.frame.paircomp <- function(x, row.names = NULL, optional = FALSE, ...,
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

## visualization of aggregated data
plot.paircomp <- function(x, off = 0.05,
  xlab = "Proportion of comparisons", ylab = "", tol.xlab = 0.05,
  abbreviate = TRUE, hue = NULL, chroma = 40, luminance = 80,
  xlim = c(0, 1), ylim = NULL, xaxs = "i", yaxs = "i", ...)
{
  ## tabulate
  tab <- summary(x)   
  ## omit NAs
  if(ncol(tab) > length(mscale(x))) tab <- tab[,-ncol(tab)]

  ## basic dimensions
  lab <- labels(x)
  nobj <- length(lab)
  ix <- which(upper.tri(diag(nobj)), arr.ind = TRUE)
  if(attr(x, "ordered")) ix <- rbind(ix, ix[,2:1])
  npc <- nrow(ix)
  has_undecided <- min(abs(mscale(x))) < 1

  ## labeling
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(lab))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  alab <- abbreviate(lab, abbreviate)  

  ## coordinates: (cumulative) proportions
  xcprob <- t(matrix(apply(tab, 1, function(z) as.vector(cumsum(z)[-length(z)]/sum(z))), ncol = nrow(tab)))
  yprob <- as.vector(rowSums(tab)/sum(tab))
  ycprob <- 1 + (npc-1) * off - c(0, cumsum(yprob[-npc] + off))

  ## HCL-based sequential palette
  if(is.null(hue)) hue <- seq(0, by = 360/nobj, length = nobj)
  hue <- rep(hue, length.out = nobj)
  chroma <- c(chroma, 0)
  luminance <- c(luminance, 95)[1:2]
  ns <- max(mscale(x)) + 1
  col <- sapply(hue, function(z) {
   rval <- seq(1, 0, length = ns)
   hcl(z, chroma[2] - diff(chroma) * rval,
     luminance[2] - diff(luminance) * rval^0.8, fixup = TRUE)
  })

  ## default arguments
  if(is.null(ylim)) ylim <- c(0, 1 + (npc-1) * off) 

  ## raw plot
  plot(0, 0, type = "n", axes = FALSE,
    xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs,
    xlab = xlab, ylab = ylab, ...)

  ## actual rectangles
  for(i in 1:npc) {
    rect(c(0, xcprob[i,]), ycprob[i] - yprob[i], c(xcprob[i,], 1), ycprob[i],
      col = c(if(has_undecided) col[,ix[i,1]] else col[-nrow(col),ix[i,1]],
      rev(col[,ix[i,2]])[-1]))
  }
  
  ## position for x axis labels
  xat <- c(xcprob[1,], 1) - diff(c(0, xcprob[1,], 1))/2
  if(any(diff(xat) < tol.xlab)) xat <- 1:ncol(tab)/(ncol(tab) + 1)
  
  ## labeling
  axis(1)
  axis(3, at = xat, labels = colnames(tab), tick = FALSE)
  axis(2, at = ycprob - yprob/2, labels = alab[ix[,1]], tick = FALSE, las = 1)
  axis(4, at = ycprob - yprob/2, labels = alab[ix[,2]], tick = FALSE, las = 1)
  box()

  ## separate within-pair orders by horizontal line
  if(attr(x, "ordered")) abline(h = ycprob[npc/2 + 1] + off/2, xpd = TRUE)

  ## return centered y coordinates invisibly
  invisible(ycprob - yprob/2)
}


## default NA handling
is.na.paircomp <- function(x) apply(is.na(as.matrix(x)), 1, all)

## FIXME: need more flexible na.* methods
