## TRAMPsamples.R (part of the TRAMPR package)

## Generator function
## TODO: fix factors here?
TRAMPsamples <- function(data, info=NULL, warn.factors=TRUE, ...) {
  if ( is.null(info) )
    info <- data.frame(sample.pk=unique(data$sample.fk))
  if ( is.null(info$species) )
    info$species <- as.character(NA)
  obj <- list(info=defactor(info, warn.factors),
              data=defactor(data, warn.factors), ...)
  class(obj) <- "TRAMPsamples"
  tidy.TRAMPsamples(obj)
}

## Identity function
is.TRAMPsamples <- function(x)
  inherits(x, "TRAMPsamples")

## Validity function
valid.TRAMPsamples <- function(x) {
  if ( !is.TRAMPsamples(x) )
    stop("Not a TRAMPsamples object")
  data.cols <- c("sample.fk", "primer", "enzyme", "size", "height")
  must.contain.cols(x$info, c("sample.pk", "species"))
  must.contain.cols(x$data, data.cols)
  if ( any(duplicated(x$info$sample.pk)) )
    stop("x$info$sample.pk must be unique")
  if ( !(is.numeric(x$info$sample.pk) && all(x$info$sample.pk > 0)) )
    stop("Numeric, positive sample.pk required")
  if ( any(is.na(x$data[data.cols])) )
    stop("NA values are not permitted for columns: ",
         paste(data.cols, collapse=", "))
  orphan <- setdiff(x$data$sample.fk, x$info$sample.pk)
  if ( length(orphan) > 0 )
    stop(sprintf("Orphaned data with no info: (sample.fk: %s)",
                 paste(orphan, collapse=", ")))
  invisible(TRUE)
}

## Other basic methods:
print.TRAMPsamples <- function(x, ...)
  cat("[TRAMPsamples object]\n")

labels.TRAMPsamples <- function(object, ...) {
  valid.TRAMPsamples(object)
  sort(object$info$sample.pk)
}

summary.TRAMPsamples <- function(object, include.info=FALSE, ...) {
  valid.TRAMPsamples(object)

  res <- object$data
  res$code <- paste(res$primer, res$enzyme, sep="_")
  res <- tapply(res$sample.fk, res[c("sample.fk", "code")], length)

  res <- res[match(object$info$sample.pk, rownames(res)),,drop=FALSE]
  rownames(res) <- object$info$sample.pk

  if ( include.info )
    cbind(object$info, res)
  else
    res
}

plot.TRAMPsamples <- function(x, sample.fk=labels(x), ...)
  for ( i in sample.fk )
    TRAMPsamples.plotone(x, i, ...)

TRAMPsamples.plotone <- function(x, sample.fk,
                                 all.samples.global=FALSE, col=1:10,
                                 xmax=NULL, mar.default=.5,
                                 mar.labels=8, cex=1) {
  sample.data <- x[sample.fk]$data

  ## Decide which enzyme/primer combinations to plot:
  d.samples <- if (all.samples.global) x$data else sample.data
  enzyme.primer <- unique(d.samples[c("enzyme", "primer")])
  enzyme.primer <-
    enzyme.primer[order(enzyme.primer$enzyme, enzyme.primer$primer),]

  ## This assertion _will_ need dealing with.
  n <- nrow(enzyme.primer)
  if ( n < 1 ) {
    warning("No data for sample: ", sample.fk)
    layout(1)
    plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
    text(.5, .5, "(No data)")
    title(main=paste("Sample:", sample.fk), outer=TRUE)
  }

  if ( length(col) < n )
    stop(sprintf("Too few colours provided (%d, %d required)",
                 length(col), n))

  sample.data$code <- classify(sample.data, enzyme.primer)

  if ( is.null(xmax) )
    xmax <- ceiling(max(x$data$size)/100)*100
  else if ( is.na(xmax) )
    xmax <- ceiling(max(sample.data$size)/100)*100

  layout(matrix(1:n, n))
  par(oma=c(4, 0, 3, 0) + mar.default, mar=c(.75, mar.labels, 0, 0),
      las=1, cex=cex)

  for ( i in seq(n) ) {
    dsub <- sample.data[sample.data$code == i,]
    plot(height ~ size, dsub, type="h", xlim=c(0, xmax),
         ylim=range(dsub$height, 0:1), xaxt="n", yaxt="n", bty="l",
         xlab="", ylab="", col=col[dsub$code])
    axis(1, labels=i == n)
    axis(2, mean(par("usr")[3:4]),
         with(enzyme.primer[i,], paste(enzyme, primer, sep="/")),
         tick=FALSE)
    if ( i == n )
      title(xlab="Fragment length", xpd=NA)
  }

  title(main=paste("Sample:", sample.fk), outer=TRUE)
}

combine.TRAMPsamples <- function(x, y, rewrite.sample.pk=FALSE, ...) {
  valid.TRAMPsamples(x)
  valid.TRAMPsamples(y)

  if ( any(labels(y) %in% labels(x)) )
    if ( rewrite.sample.pk ) {
      y$info$sample.pk <- y$info$sample.pk + max(labels(x))
      y$data$sample.fk <- y$data$sample.fk + max(labels(x))
    } else {
      stop("sample.pk conflict - see ?combine.TRAMPsamples")
    }

  x$info <- rbind2(x$info, y$info)
  x$data <- rbind2(x$data, y$data)

  extra <- setdiff(names(y), c("info", "data"))
  if ( length(extra) > 0 )
    warning("Additional objects in 'y' were ignored: ",
            paste(dQuote(extra), collapse=", "))

  tidy.TRAMPsamples(x)
}

tidy.TRAMPsamples <- function(x) {
  valid.TRAMPsamples(x)
  x$info <- x$info[order(x$info$sample.pk),]
  x$data <- x$data[do.call(order, x$data),]
  rownames(x$info) <- seq(length=nrow(x$info))
  rownames(x$data) <- seq(length=nrow(x$data))
  x
}

## _Basic_ indexing method.
"[.TRAMPsamples" <- function(x, i, na.interp=TRUE, ...) {
  if ( is.numeric(i) ) {
    if ( !all(i %in% labels(x)) )
      stop("Unknown samples: ",
           paste(i[!(i %in% labels(x))], collapse=", "))
  } else if ( is.logical(i) ) {
    if ( length(i) != nrow(x$info) )
      stop("Logical index of incorrect length")
    i[is.na(i)] <- na.interp
    i <- x$info$sample.pk[i]
  }  else stop("Invalid index type")

  ## These generate NOTEs in R CMD CHECK
  x$info <- subset(x$info, x$info$sample.pk %in% i)
  x$data <- subset(x$data, x$data$sample.fk %in% i)
  ##x$info <- x$info[x$info$knowns.pk %in% i,]
  ##x$data <- x$data[x$data$knowns.pk %in% i,]
  
  tidy.TRAMPsamples(x)
}
