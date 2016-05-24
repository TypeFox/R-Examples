## TRAMP.R (part of the TRAMPR package)

## Using rowSums() is about 200 times faster than apply(, 1:2, sum,
## ...), but there is no rowMax() function.

## When using max(x, na.rm=TRUE), if all x are NA, then you get a
## warning and -Inf as the value.  Max.na.rm() (within TRAMP) returns
## NA in this case.
TRAMP <- function(samples, knowns, accept.error=1.5, min.comb=4,
                  method="maximum") {
  diffsmatrix <- create.diffsmatrix(samples, knowns)
  if ( is.na(min.comb) )
    min.comb <- dim(diffsmatrix)[3]

  n <- rowSums(!is.na(diffsmatrix), dims=2)

  Max.na.rm <- function(x)
    if ( all(is.na(x)) ) NA else max(x, na.rm=TRUE)

  method <- match.arg(method, c("maximum", "euclidian", "manhattan"))
  if ( method == "maximum" )
    error <- apply(abs(diffsmatrix), 1:2, Max.na.rm)
  else if ( method == "euclidian" )
    error <- sqrt(rowSums(diffsmatrix^2, TRUE, 2))/n
  else if ( method == "manhattan" )
    error <- rowSums(abs(diffsmatrix), TRUE, 2)/n

  presence.ign <-
    data.frame(sample.fk=numeric(0), knowns.fk=numeric(0))

  res <- list(presence      = n >= min.comb & error < accept.error,
              presence.ign  = presence.ign,
              error         = error,
              n             = n,
              diffsmatrix   = diffsmatrix,
              enzyme.primer = attr(diffsmatrix, "enzyme.primer"),
              samples       = samples,
              knowns        = knowns,
              accept.error  = accept.error,
              min.comb      = min.comb,
              method        = method)
  class(res) <- "TRAMP"
  res
}

rebuild.TRAMP <- function(x) {
  if ( !is.TRAMP(x) )
    stop("x must be a TRAMP object")
  TRAMP(x$samples, x$knowns, accept.error=x$accept.error,
        min.comb=x$min.comb, method=x$method)
}

## Create the main "diffmatrix"; a 3d matrix with the distances
## between peaks in the knowns database and the sample data.

## For each sample, this computes the maximum (absolute) distance
## between

## While the minimim absolute distance is calculated, the sign of the
## distance is preserved (negative where 

## For ~450 samples and ~300 knowns, this takes ~20s on  3GHz
## machine.  However, it need be run only once.

## In detail;

## TODO: This documentation is painfully old, and needs to be updated
## in places to reflect new data structures.

## Determine the number and indices of the unique samples in the
## rundata (samples/n.sample) and the knowns database
## (knowns.sample.fk/n.known) [determined by `$sample.fk'].

## Collect all unique enzyme/primer combinations used in the rundata
## (`enzyme.primer'); each row is a different combination, and there
## are `n.comb' combinations.

## In the rundata, create a column `$code', indicating which row in
## enzmye.primer the combination corresponds to (this is just to
## simplify a later step).

## Split the knowns data into a list of n.comb elements, each
## corresponding to peak sizes for the different primer/enzyme
## combinations.  Each element is of length n.sample, and in the order
## defined by knowns.sample.fk; where a known lacks a peak for a
## particular enzyme/primer combination, an NA is used, so that all
## elements in `knowns.split' have n.knowns elements, in the order of
## knowns$info$sample.fk

## The main loop works over each sample, then each enzyme/primer
## combination.
## For each it calculates the distance between every peak in the data
## (for this sample+enzmye+primer) and every peak in the
## knowns database (for this enzyme+primer), and takes the minimum
## across different sample peaks:
##    ** this gives the _smallest_ distance between a known peak and
##       _any_ sample peak.
## This vector is of length n.knowns.  While the smallest absolute
## distance is taken, the sign is preserved (negative where the
## closest sample peak was less than the known peak, positive where
## greater).

## For each individual this generates an n.known x n.comb matrix,
## across the n.comb different primer/enzyme combinations.  It is
## against this matrix that identification methods should work; using
## the sum of the differences, the sum of squares, the largest
## absolute difference, etc.  For a known to be matched, all peaks
## shared by the sample and the knowns database should be "close
## enough".  The decision of "close enough" may not be clear,
## however.

## Across all samples, this generates an n.sample x n.known x n.comb
## array, where
##   m[i,,]
## gives the 2-way matrix described above for the ith individual.

## To record which of the n.comb "rows" correspond to which
## primer/enzyme combination, an attribute "primer.enzyme" is included
## with the output array; it is a data.frame of two columns (primer
## and enzyme), such that the enzyme and primer used in
##    m[,,i]
## is
##    enzyme.primer$enzyme[i]
##    enzyme.primer$primer[i]
## respectively.

## OK, this is only working where enzyme/primer combinations are
## shared across the data and the knowns.

## Previously I've been looping through each sample, then
## enzyme/primer.  Going through enzyme/primers first might be
## better:
create.diffsmatrix <- function(samples, knowns) {
  valid.TRAMPsamples(samples)
  valid.TRAMPknowns(knowns)

  sample.fk <- labels(samples)
  knowns.fk <- labels(knowns)
  n.samples <- length(sample.fk)
  n.knowns  <- length(knowns.fk)

  samples.data <-
    samples$data[c("sample.fk", "enzyme", "primer", "size")]
  knowns.data <-
    knowns$data[c("knowns.fk", "enzyme", "primer", "size")]

  enzyme.primer <- unique(knowns.data[c("primer", "enzyme")])
  n.comb <- nrow(enzyme.primer)
  if ( !(n.samples > 0 && n.knowns > 0) )
    stop("At least one sample and known is required")
  if ( is.null(n.comb) || n.comb < 1 )
    stop("At least one enzyme/primer combination required")

  knowns.data$code <- classify(knowns.data, enzyme.primer)
  samples.data$code <- classify(samples.data, enzyme.primer)
  samples.data <- samples.data[!is.na(samples.data$code),]

  res <- vector("list", n.comb)
  for ( comb in 1:n.comb ) {
    s.sub <- samples.data[samples.data$code == comb,]
    if ( nrow(s.sub) < 1 ) {
      res[[comb]] <- matrix(NA, n.knowns, n.samples)
      next
    }

    k.sub <- knowns.data[knowns.data$code == comb,]
    k <- k.sub$size
    m <- matrix(NA, length(k), n.samples)
    for ( sample in 1:n.samples ) {
      s <- s.sub$size[s.sub$sample.fk == sample.fk[sample]]
      if ( length(s) > 0 )
        m[,sample] <- apply(outer(s, k, "-"), 2, absolute.min)
    }

    res[[comb]] <- m[match(knowns.fk, k.sub$knowns.fk),]
  }

  if ( nrow(unique(t(sapply(res, dim)))) != 1 )
    stop("Failed to construct diffsmatrix")

  m <- array(unlist(res), c(n.knowns, n.samples, n.comb),
             dimnames=list(known=knowns.fk, sample=sample.fk,
               enzyme.primer=NULL))
  m <- aperm(m, c(2,1,3))

  rownames(enzyme.primer) <- 1:n.comb
  attr(m, "enzyme.primer") <- enzyme.primer
  m
}

is.TRAMP <- function(x)
  inherits(x, "TRAMP")

print.TRAMP <- function(x, ...)
  cat("[TRAMP object]\n")

## I'd rather use 'x' here as the first argument, but summary
## dispatches on 'object'.
summary.TRAMP <- function(object, name=FALSE, grouped=FALSE,
                          ignore=FALSE, ...) {
  m <- object$presence
  knowns <- object$knowns
 
  if ( ignore )
    m[cbind(match(object$presence.ign$sample.fk, rownames(m)),
            match(object$presence.ign$knowns.fk, colnames(m)))] <-
              FALSE

  if ( grouped ) {
    m <- t(apply(m, 1, tapply, knowns$info$group.strict, any))
    if ( name )
      colnames(m) <- 
        knowns$info$group.name[match(colnames(m),
                                     knowns$info$group.strict)]
  } else if ( name )
    colnames(m) <- 
      knowns$info$species[match(colnames(m), labels(object$knowns))]
  m
}

remove.TRAMP.match <- function(x, sample.fk, knowns.fk) {
  if ( !is.TRAMP(x) )
    stop(dQuote(x), " must be a TRAMP object")
  if ( length(sample.fk) != 1 || length(knowns.fk) != 1 )
    stop("sample.fk and knowns.fk must be of length 1")
  if ( !sample.fk %in% labels(x$samples) )
    stop("Unknown sample.fk ", dQuote(sample.fk))
  if ( !knowns.fk %in% labels(x$knowns) )
    stop("Unknown knowns.fk ", dQuote(knowns.fk))
  if ( any(x$presence.ign$sample.fk == sample.fk &
           x$presence.ign$knowns.fk == knowns.fk) )
    stop("sample.fk/knowns.fk combination already ignored")
  m <- x$presence
  if ( !m[as.character(sample.fk), as.character(knowns.fk)] )
    stop(sprintf("%d/%d (sample.fk/knowns.fk) is not a match",
                 sample.fk, knowns.fk))
  x$presence.ign[nrow(x$presence.ign)+1,] <-
    list(sample.fk, knowns.fk)
  x
}

## Plotting method:
plot.TRAMP <- function(x, sample.fk=labels(x$samples), ...)
  for ( i in sample.fk )
    TRAMP.plotone(x, i, ...)

## Plot a single TRAMP fit.

## Layout model:
## The outside of the plot is padded with 'mar.default' lines of space
## Three lines each is reserved at the top of the device (outer
##   margin) for the title.
## Three lines is reserved below each subplot for the x-axis
## The proportion 'p.labels' of the available horizontal space is
##   reserved for labels.

## Knowns that match this data, grouped appropriately if grouped
## is TRUE.  'ylab' is the vector of knowns that we are plotting
## against.

## The issue here is how to deal with differences in presence of
## enzymes and primers between the samples and knowns.

## Introduce two new arguments:
##   all.samples all.knowns
##   FALSE       FALSE      not allowed
##   TRUE        FALSE      sample combinations
##   FALSE       TRUE       knowns combinations
##   TRUE        TRUE       union

## FALSE/FALSE could compute the intersection, but I can't be bothered
## implementing that.

## However, for all.samples, should be scope be global (so that the
## same set of plots is used for all cases) or local (so that only
## cases present in the _current_ sample are plotted).  This is
## decided by the value of 'all.sample.global'

TRAMP.plotone <- function(x, sample.fk, grouped=FALSE, ignore=FALSE,
                          all.knowns=TRUE, all.samples=FALSE,
                          all.samples.global=FALSE, col=1:10,
                          pch=if (grouped) 15 else 16, xmax=NULL,
                          horiz.lines=TRUE,  mar.default=.5, p.top=.5,
                          p.labels=1/3, cex.axis=NULL,
                          cex.axis.max=1) {
  if ( !sample.fk %in% labels(x$samples) )
    stop("No such sample.fk ", dQuote(sample.fk))

  op <- par(no.readonly=TRUE)
  on.exit(par(op))

  sample.data <- x$samples[sample.fk]$data
  
  ## Construct the enzyme/primers, and order by enzyme, then primer.
  if ( !all.samples && !all.knowns )
    stop("One of all.samples and all.knowns must be TRUE")
  if ( all.samples ) {
    d.samples <- if (all.samples.global)
      x$samples$data else sample.data
    enzyme.primer <- unique(d.samples[c("enzyme", "primer")])
  }
  if ( all.knowns ) {
    if ( all.samples )
      enzyme.primer <-
        unique(rbind(x$enzyme.primer, enzyme.primer))
    else
      enzyme.primer <- x$enzyme.primer
  }
  enzyme.primer <-
    enzyme.primer[order(enzyme.primer$enzyme, enzyme.primer$primer),]

  n <- nrow(enzyme.primer)
  stopifnot(n > 0)

  if ( length(col) < n )
    stop(sprintf("Too few colours provided (%d, %d required)",
                 length(col), n))

  sample.data$code <- classify(sample.data, enzyme.primer)
  sample.data <- sample.data[!is.na(sample.data$code),]

  match.fk <- 
    as.integer(names(which(summary(x, ignore=ignore)[
                                        as.character(sample.fk),])))
  matches.info <- subset(x$knowns$info,
                         x$knowns$info$knowns.pk %in% match.fk)

  if ( !nrow(matches.info) )
    ylab <- "(No matches)"
  else {
    i <- is.na(matches.info$species)
    if ( any(i) )
      matches.info$species[i] <-
        sprintf("[Unknown sample: %d]", matches.info$knowns.pk[i])
    
    if ( grouped ) {
      matches.info$group.name <- factor(matches.info$group.name)
      matches.info <- matches.info[order(matches.info$group.name),]
      ylab <- levels(matches.info$group.name)
    } else {
      matches.info <- matches.info[order(matches.info$species),]
      ylab <- matches.info$species
    }
  }

  matches.data <-
    subset(x$knowns$data,
           x$knowns$data$knowns.fk %in% matches.info$knowns.pk)
  matches.data$code <- classify(matches.data, enzyme.primer)

  matches.data$row <-
    match(matches.data$knowns.fk, matches.info$knowns.pk)
  if ( grouped )
    matches.data$row <-
      as.integer(matches.info$group.name)[matches.data$row]

  ## If xmax is not given, take the maximum observed value across the
  ## dataset, rounded to the nearest hundred (rounded up only).
  if ( is.null(xmax) )
    xmax <- ceiling(max(x$samples$data$size)/100)*100
  else if ( is.na(xmax) )
    xmax <- ceiling(max(sample.data$size)/100)*100

  ## Start plotting; set aside p.top of the vertical space for the
  ## knowns vs. data plot, with the rest for peaks by enzyme/primer
  ## combination.  Set aside p.labels of the horizontal space on the
  ## LHS for group labels (i.e. ylab).
  layout(matrix(1:(n+1), n+1), heights=c(p.top, rep((1-p.top)/n, n)))
  par(oma=c(4, 0, 3, 2) + mar.default, mar=c(.75, 0, 0, 0))
  mai <- par("mai")
  mai[2] <- par("fin")[1]*p.labels
  par(mai=mai)

  plot(NA, xlim=c(0, xmax), ylim=c(length(ylab), 0)+0.5, type="n",
       bty="l", xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  if ( horiz.lines && nrow(matches.info) )
    abline(h=1:length(ylab), col="gray", lty=3)
  abline(v=sample.data$size, lty=2, col=col[sample.data$code])

  ## cex.axis must be calculated _after_ the plot area is set up.
  if ( is.null(cex.axis) )
    cex.axis <- fit.cex.yaxis(ylab, cex.axis.max)
  par(cex.axis=cex.axis, tcl=-cex.axis/2, mgp=c(3, cex.axis, 0),
      las=1)
  cex.lab <- if ( cex.axis > 1 ) cex.axis else 1
  
  axis(1, labels=FALSE)
  axis(2, 1:length(ylab), ylab, tick=nrow(matches.info))

  points(matches.data$size, matches.data$row, pch=pch,
         cex=1.5, col=col[matches.data$code])

  ## Display peaks for each enzyme.primer combination.
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
      title(xlab="Fragment length", xpd=NA, cex.lab=cex.lab)
  }

  title(main=paste("Sample:", sample.fk), outer=TRUE,
        cex.main=cex.lab)
  mtext(c("Matches", "Enzyme/Primer combinations"), 4, 1,
        at=c(.75, .25), outer=TRUE, las=0, cex=.8*cex.lab)

  ## Not sure if we actually want to trust this, since it could be
  ## subset by codes...
  invisible(list(idx=sample, sample.data=sample.data,
                 matches.info=matches.info,
                 matches.data=matches.data))
}

fit.cex.yaxis <- function(ylab, cex.axis.max) {
  cheight <- par("csi")
  cex.v <- .9/cheight
  cex.h <- par("mai")[2]/(1.5*cheight + max(strwidth(ylab, "inches")))
  min(cex.v, cex.h, cex.axis.max)
}

## New update function (complete rewrite)

## TODO: Because I need to trap errors nicely, I'm should wrap the
## actual bits that do anything (ADD, IGNORE and SAVE) in functions
## (since R has no block-exception syntax), and then do try() on them.
## However, without a call/cc variant it's fairly hard to control how
## that works?

## There are possibilities of failing in all of these:
## SAVE: (e.g. if saving fails due to an open file, full file system,
##   etc.)
## ADD: Add knowns already warns nicely, but it's still possible this
##   could fail.
## IGNORE: Really shouldn't fail, but might.

## Grr.  Let's assume no errors occur, to keep things simple:
update.TRAMP <- function(object, sample.fk=labels(object$samples),
                         grouped=FALSE, ignore=TRUE,
                         delay.rebuild=FALSE,
                         default.species=NULL,
                         filename.fmt="TRAMP_%d.pdf", ...) {
  ## grouped=TRUE is not allowed, because it interferes in a not-nice
  ## way with ignoring matches.
  if ( grouped ) .NotYetUsed("grouped != FALSE")
  if ( !ignore ) .NotYetUsed("grouped != TRUE")

  if ( getOption("warn") == 0 ) {
    oo <- options(warn=1)
    on.exit(options(oo))
  }

  known.added <- FALSE
  rebuild <- !delay.rebuild

  choices <- c("Add a sample to the knowns database",
               "Mark a match to be ignored",
               "Save the current plot as a PDF",
               "Move to the next plot",
               "Quit")
  ADD <- 1; IGNORE <- 2; SAVE <- 3; NEXT <- 4; QUIT <- c(0, 5)

  for ( i in sample.fk ) {
    repeat {
      plot.obj <- plot(object, i, grouped=grouped, ignore=ignore, ...)
      action <- menu(choices)
    
      if ( action == ADD )
        object <- add.known(object, i, rebuild,
                            default.species=default.species)
      else if ( action == IGNORE ) {
        matches <- plot.obj$matches.info[c("knowns.pk", "species")]
        if ( nrow(matches) > 0 ) {
          ign <- menu(matches$species,
                      title="Select match to exclude (0 cancels)")
          object <- remove.TRAMP.match(object, i,
                                       matches$knowns.pk[ign])
        } else
        warning("No matches to ignore!")
      } else if ( action == SAVE ) {
        filename <- sprintf(filename.fmt, i)
        cat(sprintf("\tPlot saved as %s\n", filename))
        dev.copy(pdf, file=filename)
        dev.off()
      } else if ( action %in% c(NEXT, QUIT) )
        break
    }
    if ( action %in% QUIT )
      break
  }

  if ( known.added && delay.rebuild ) {
    cat("Rebuilding TRAMP object\n")
    object <- rebuild.TRAMP(object)
  }
  
  object
}
