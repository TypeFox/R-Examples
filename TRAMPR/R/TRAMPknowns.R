## TRAMPknowns.R (part of the TRAMPR package)

## Generator function
TRAMPknowns <- function(data, info, cluster.pars=list(),
                        file.pat=NULL, warn.factors=TRUE, ...) {
  cluster.pars.defaults <- list(dist.method="maximum",
                                hclust.method="complete",
                                cut.height=2.5)
  ok <- names(cluster.pars) %in% names(cluster.pars.defaults)
  if ( any(!ok) )
    stop("Invalid cluster.pars entries given: ",
         paste(names(cluster.pars)[!ok], collapse=", "))
  cluster.pars.defaults[names(cluster.pars)] <- cluster.pars

  obj <- list(info=defactor(info, warn.factors),
              data=defactor(data, warn.factors),
              cluster.pars=cluster.pars.defaults,
              file.pat=file.pat, ...)
  class(obj) <- "TRAMPknowns"
  tidy.TRAMPknowns(obj)
}

## Identity function
is.TRAMPknowns <- function(x)
  inherits(x, "TRAMPknowns")

## Validity function
## TODO: fix factors here (or even test for presence)?
## Classes:
##   info:
##     knowns.pk: numeric (tested)
##     species: character/NA
##   data:
##     knowns.fk: numeric
##     primer: character
##     enzyme: character
##     size: numeric
valid.TRAMPknowns <- function(x) {
  if ( !is.TRAMPknowns(x) )
    stop("Not a TRAMPknowns object")

  data.cols <- c("knowns.fk", "primer", "enzyme", "size")
  must.contain.cols(x$info, c("knowns.pk", "species"))
  must.contain.cols(x$data, data.cols)
  if ( any(duplicated(x$info$knowns.pk)) )
    stop("x$info$knowns.pk must be unique")
  if ( !(is.numeric(x$info$knowns.pk) && all(x$info$knowns.pk > 0)) )
    stop("Numeric, positive knowns.pk required")
  if ( any(is.na(x$data[data.cols])) )
    stop("NA values are not permitted for columns: ",
         paste(data.cols, collapse=", "))

  data.info <- x$data[c("knowns.fk", "primer", "enzyme")]
  if ( any(duplicated(data.info)) ) {
    i <- duplicated(data.info)
    bad <- x$data[!is.na(classify(x$data, data.info[i,])),]
    bad <- bad[do.call(order, bad),]
    rownames(bad) <- paste("  error:", rownames(bad), sep="")
    on.exit(print(bad))
    stop("Duplicate enzyme/primer combinations within a knowns.fk:")
  }

  ## Check for orphaned data and data-less knowns.  Both are invalid
  ## (i.e. a known _must_ have a peak profile, otherwise it is not a
  ## known.  A peak _must_ be associated with a known, otherwise it is
  ## junk data.
  no.data <- setdiff(x$info$knowns.pk, x$data$knowns.fk)
  orphan <- setdiff(x$data$knowns.fk, x$info$knowns.pk)

  if ( length(no.data) > 0 )
    stop(sprintf("Knowns in database with no data: (knowns.pk: %s)",
                 paste(no.data, collapse=", ")))
  if ( length(orphan) > 0 )
    stop(sprintf("Orphaned data with no info: (knowns.fk: %s)",
                 paste(orphan, collapse=", ")))
    
  invisible(TRUE)
}

## This performs some basic cleanup after directly manipulating a
## knowns database (i.e. a TRAMPknowns object)
## This is run after adding knowns (combine.TRAMPknowns() and
## add.known()), and for the initial creation of a TRAMPknowns object
## (TRAMPknowns()).

## The function just sorts the data.frames, keeps the row.names
## in order (may change for R 2.4.0) and reclusters the data.
tidy.TRAMPknowns <- function(x) {
  valid.TRAMPknowns(x)
  x$info <- x$info[order(x$info$knowns.pk),]
  x$data <- x$data[do.call(order, x$data),]
  rownames(x$info) <- seq(length=nrow(x$info))
  rownames(x$data) <- seq(length=nrow(x$data))
  group.knowns(x)
}

## Grouping functions:
## This is the only function that will be exported.
group.knowns <- function(x, ...)
  UseMethod("group.knowns")
  
group.knowns.TRAMPknowns <- function(x, dist.method, hclust.method,
                                     cut.height, ...) {
  valid.TRAMPknowns(x)
  name.knowns(collapse.knowns(cluster.knowns(x, dist.method,
                                             hclust.method,
                                             cut.height)))
}

group.knowns.TRAMP <- function(x, ...) {
  x$knowns <- group.knowns(x$knowns, ...)
  x
}

## Clustering function:
## (the missing() tells if the argument is missing in the calling
## function (group.knowns()).
cluster.knowns <- function(x, dist.method, hclust.method,
                           cut.height) {
  pars <- x$cluster.pars
  if ( !missing(cut.height) ) pars$cut.height <- cut.height
  if ( !missing(dist.method) ) pars$dist.method <- dist.method
  if ( !missing(hclust.method) ) pars$hclust.method <- hclust.method
  x$cluster.pars <- pars

  peaks <- summary(x)
  if ( nrow(peaks) > 1 ) {
    cluster <- hclust(dist(peaks, pars$dist.method), pars$hclust.method)
    groups <- cutree(cluster, h=pars$cut.height)
    x$info$group.cluster <-
      groups[match(x$info$knowns.pk, names(groups))]
    x$cluster <- cluster
  } else if ( nrow(peaks) == 1 ) {
    x$cluster <- NA
    x$info$group.cluster <- x$info$group.strict <- 1
    x$info$group.name <- x$info$species
  } else stop("Cannot cluster empty TRAMPknowns object")
  x
}

## A function to collapse both species and group.clusters, to create
## "strict" groups.

## This is iterative, and should take many fewer iterations than there
## are different knowns types.  For each iteration, collapse over both
## the species and knowns.group.fk indices.
## The function 'f' does the collapsing at each iteration.

## If 'species' is NA (an unidentified pattern), then it will be given
## its knowns.pk value (after checking that no species has that
## value).  This is a bit of a kludge, however.
collapse.knowns <- function(x) {
  info <- x$info
  must.contain.cols(info, "group.cluster")

  species <- info$species
  i <- is.na(species)
  if ( any(info$knowns.pk[i] %in% species[!i]) )
    stop("knowns.pk values present in species!")
  species[i] <- info$knowns.pk[i]

  f <- function(x, index) {
    index <- factor(index)
    as.vector(tapply(x, index, min)[as.integer(index)])
  }

  groups <- 1:nrow(info)
  groups.old <- integer(0)
  for ( n in 1:nrow(info) ) {
    groups <- f(f(groups, species), info$group.cluster)
    if ( identical(groups, groups.old) )
      break
    else
      groups.old <- groups
  }

  info$group.strict <- as.integer(factor(groups))
  x$info <- info
  x
}

name.knowns <- function(x) {
  info <- x$info
  must.contain.cols(info, "group.strict")

  ## Given 'i' as a vector of row-indices in info:
  ## (Take only the first 'n' unknown groups if there are more than n
  ## available.  Always take all species).
  generate.name <- function(i) {
    if ( any(!is.na(info$species[i])) )
      paste(sort(unique(na.omit(info$species[i]))), collapse=", ")
    else
      sprintf("Unknown cluster (knowns.pk: %s)",
              paste(sort(info$knowns.pk[i]), collapse=", "))
  }
  info$group.name <- 
    as.vector(tapply(1:nrow(info), info$group.strict,
                     generate.name)[info$group.strict])
  x$info <- info
  x
}

## Other basic methods:
print.TRAMPknowns <- function(x, ...)
  cat("[TRAMPknowns object]\n")

labels.TRAMPknowns <- function(object, ...) {
  valid.TRAMPknowns(object)
  sort(object$info$knowns.pk)
}

summary.TRAMPknowns <- function(object, include.info=FALSE, ...) {
  valid.TRAMPknowns(object)
  res <- object$data
  res$code <- paste(res$primer, res$enzyme, sep="_")
  res <- reshape(res[c("knowns.fk", "code", "size")],
                 direction="wide", idvar="knowns.fk",
                 timevar="code")
  colnames(res) <- sub("^size\\.", "", colnames(res))

  res <- res[match(object$info$knowns.pk, res$knowns.fk),-1]
  rownames(res) <- object$info$knowns.pk

  if ( include.info )
    cbind(object$info, res)
  else
    res
}


## Plotting function for the knowns.

## The idea here is to display the clustering and profiles of the
## knowns, to inspect how things are being joined together.

## The plot will be fairly unreadable on the screen, but should print
## out OK.
plot.TRAMPknowns <- function(x, cex=1, name="species", pch=1,
                             peaks.col, p=.02, group.clusters=TRUE,
                             groups.col=1:4, grid.by=5,
                             grid.col="gray", widths=c(1, 2, 1),
                             ...) {
  valid.TRAMPknowns(x)
  if ( length(labels(x)) < 2 )
    stop("Cannot draw a plot for less than two knowns. Sorry.")
  cluster <- as.dendrogram(x$cluster)
  peaks <- summary(x)

  cluster.order <- match(labels(cluster), rownames(peaks))
  n <- nrow(peaks)
  if ( missing(peaks.col) )
    peaks.col <- 1:ncol(peaks)
  else if ( length(peaks.col) < ncol(peaks) )
    stop(sprintf("Too few colours provided (%d, %d required)",
                 length(peaks.col), n))
  group.cluster <- x$info$group.cluster[cluster.order]
  group.strict <- x$info$group.strict[cluster.order]

  grid.at <- if ( !is.na(grid.by) && !is.na(grid.col) &&
                 ceiling(grid.by/2) < n )
    seq(ceiling(grid.by/2), n, by=grid.by) else grid.at <- NA

  ## Build strings of species or group names:
  name <- match.arg(name, c("species", "group.name"))
  if ( name == "species" ) {
    i <- is.na(x$info$species)
    if ( any(i) )
      x$info$species[i] <-
        sprintf("Unknown sample %d", x$info$knowns.pk[i])
  }
  f <- if ( name == "species" )
    function(x) paste(sort(unique(x)), collapse=", ") else unique
  lab <- tapply(x$info[[name]][cluster.order], group.cluster, f)
  stopifnot(all(sapply(lab, length) == 1))

  layout(matrix(c(1,3,2), 1), widths=widths)
  par(mar=c(5.1, 0, 2, 0), oma=rep(1, 4))

  ## (1) Dendrogram:
  plot(cluster, leaflab="none", horiz=TRUE, xlab="Distance",
       edgePar=list(lwd=cex))
  mtext("Dendrogram", 3, line=1)
  usr <- par("usr")

  ## (2) Peak profiles:
  plot(NA, type="n", bty="n", xlab="Size (bp)", yaxt="n",
       ylab="", xlim=c(0, ceiling(max(peaks, na.rm=TRUE)/100)*100),
       ylim=usr[3:4], yaxs="i")
  abline(h=grid.at, col=grid.col, lwd=cex/2)
  matpoints(peaks[cluster.order,], 1:n, col=peaks.col[1:ncol(peaks)],
            pch=pch, cex=cex, lwd=cex)
  mtext("Peak profiles", 3, line=1)

  ## (3) Assorted information:
  plot(NA, type="n", xlim=0:1, ylim=par("usr")[3:4], xaxs="i",
       yaxs="i", axes=FALSE, xlab="", ylab="")
  abline(h=grid.at, col=grid.col, lwd=cex/2)

  w <- c(max(strwidth(labels(cluster), "user", cex)),
         max(strwidth(group.cluster, "user", cex)),
         max(strwidth(group.strict, "user", cex)))

  ## (a) knowns.pk, group.cluster and group.strict labels.
  text(p,                 1:n, labels(cluster), adj=0, cex=cex)
  text(p*2 + w[1],        1:n, group.cluster,   adj=0, cex=cex)
  text(p*3 + sum(w[1:3]), 1:n, group.strict,    adj=1, cex=cex)

  ## (b) group.cluster brackets:
  same.grp <- which(diff(group.cluster) == 0)
  b.x0 <- p*3.5 + sum(w[1:3])
  b.x1 <- p*4   + sum(w[1:3])
  b.y0 <- tapply(same.grp, group.cluster[same.grp], min) - 1/3
  b.y1 <- tapply(same.grp, group.cluster[same.grp], max) + 1/3 + 1
  if ( length(same.grp) > 0 ) {
    segments(b.x1, b.y0, b.x1, b.y1, lwd=cex)
    segments(b.x1, b.y0, b.x0, b.y0, lwd=cex)
    segments(b.x1, b.y1, b.x0, b.y1, lwd=cex)
  }

  ## (c) group.strict joins (j=join)
  ## Calculate the "middle" of each group (used by both joins and
  ## labels).
  group.y <- tapply(1:n, group.cluster, function(x) sum(range(x))/2)

  if ( group.clusters ) {
    joint <- tapply(group.cluster, group.strict, unique)
    joint <- joint[sapply(joint, length) > 1]
    j.x <- ifelse(tapply(1:n, group.cluster, length) > 1, b.x1, b.x0)
    j.x <- lapply(joint, function(x) j.x[as.character(x)])
    j.y <- lapply(joint, function(x) group.y[as.character(x)])
    j.x1 <- b.x1
    groups.col <- rep(groups.col, length=length(joint))

    for ( j in order(sapply(j.y, "[", 1)) ) {
      j.x1 <- j.x1 + p*2/3
      segments(j.x1, min(j.y[[j]]), j.x1, max(j.y[[j]]), lwd=cex,
               col=groups.col[j])
      segments(j.x1, j.y[[j]], j.x[[j]], j.y[[j]], lwd=cex,
               col=groups.col[j])
    }
  } else j.x1 <- b.x1

  ## (d) species/group.name labels
  lab.x <- j.x1 + p
  wid.total <- 1 - p - lab.x
  if ( wid.total < 0 )
    stop("Not enough space for species names!")

  wid <- strwidth(lab, "user", cex)

  ## Determine the number of rows available for each label:  Multi-row
  ## labels are only available to clustered groups, but these will
  ## generally have the longest names.  The exception is where a
  ## clustered group joins a named group.  Wrap the labels using my
  ## function labwrap() (based on strwrap, in util.R), and for any
  ## labels that are still too long, truncate with "...").
  lab.wrap <- labwrap(lab, wid.total, exdent=2, cex=cex)
  nlines <- pmax(floor(tapply(1:n, group.cluster, length)/
                       (1.3 * strheight(strheight(lab, cex=cex)))), 1)
  for ( i in which(sapply(lab.wrap, length) > nlines) ) {
    tmp <- lab.wrap[[i]][seq_len(nlines[i])]
    tmp[nlines[i]] <- paste(tmp[nlines[i]], "...", sep="")
    lab.wrap[[i]] <- tmp
  }
  lab.wrap <- sapply(lab.wrap, paste, collapse="\n")
  text(lab.x, group.y, lab.wrap, cex=cex, adj=0)

  ## (e) decorate axis
  at.axis <- c(p + w[1]/2, 2*p + w[1] + w[2]/2,
               3*p + sum(w[1:2]) + w[3]/2, lab.x + 3*p)
  axis(1, at.axis,
       c("knowns.pk", "group.cluster", "group.strict", name),
       mgp=c(0,0,0), tcl=0, las=2, tick=FALSE)

  invisible()
}

add.known <- function(x, ...)
  UseMethod("add.known")

add.known.TRAMP <- function(x, sample.fk, rebuild=TRUE, ...) {
  x$knowns <- add.known(x$knowns, x$samples, sample.fk, ...)
  if ( rebuild )
    x <- rebuild.TRAMP(x)
  x
}

add.known.TRAMPknowns <- function(x, samples, sample.fk, prompt=TRUE,
                                  default.species=NULL, ...) {
  ## valid.TRAMPknowns() and valid.TRAMPsamples() will be called by
  ## labels(x) and labels(samples), respectively.
  if ( sample.fk %in% labels(x) ) {
    warning("This sample already in database (not added)")
    return(x)
  }

  if ( !(sample.fk %in% labels(samples)) )
    stop("No such sample.fk ", dQuote(sample.fk))

  if ( is.null(default.species) )
    default.species <-
      samples$info$species[samples$info$sample.pk == sample.fk]

  if ( interactive() && prompt ) {
    cat("Enter a species name for this new known, ",
        "or press <enter> for ", sQuote(default.species), "\n",
        sep="")
    new.name <- read.string("[New name]> ", default.species)
  } else new.name <- default.species

  new.data <-
    find.max.peaks(samples$data[samples$data$sample.fk == sample.fk,])
  if ( !all(is.na(new.data$ratio)) )
    warning("Mutiple peaks per enzyme/primer combination for ",
            "sample ", sample.fk, ", taking largest peaks only")
  names(new.data)[names(new.data) == "sample.fk"] <- "knowns.fk"
  new.data <- new.data[c("knowns.fk", "primer", "enzyme", "size")]

  x$info[nrow(x$info) + 1, c("knowns.pk", "species")] <-
    list(sample.fk, new.name)
  x$data <- rbind2(x$data, new.data)
  x <- tidy.TRAMPknowns(x)
  write.TRAMPknowns(x, warn=FALSE)
  x
}


combine <- function(x, y, ...)
  UseMethod("combine")

combine.TRAMP <- function(x, y, rebuild=TRUE, ...) {
  x$knowns <- combine(x$knowns, y, ...)
  if ( rebuild )
    x <- rebuild.TRAMP(x)
  x
}

combine.TRAMPknowns <- function(x, y, rewrite.knowns.pk=FALSE,
                                ...) {
  ## valid.TRAMPknowns() will be called by labels(x) and labels(y).

  ## TODO: This is very crude, and should be improved.
  if ( any(labels(y) %in% labels(x)) )
    if ( rewrite.knowns.pk ) {
      y$info$knowns.pk <- y$info$knowns.pk + max(labels(x))
      y$data$knowns.fk <- y$data$knowns.fk + max(labels(x))
    } else {
      stop("knowns.pk conflict - see ?combine.TRAMPknowns")
    }

  x$info <- rbind2(x$info, y$info)
  x$data <- rbind2(x$data, y$data)

  x <- tidy.TRAMPknowns(x)
  write.TRAMPknowns(x, warn=FALSE)
  x
}

## Below here unchecked, really.

## Collect the sample.fk/primer/enzyme combinations in the same
## order that peak sizes/heights will be generated.
## This just grabs the first element of sample.fk, primer and
## enzyme, since these are known not to change within levels of
## 'i'.


## The inner function 'f' does most of the work: Given the indices
## (within ) of data for one enzyme/primer combination, calculate the
## location (size) and heights of the two largest peaks.

## Only count cases where height1/height2 > min.ratio.  If height2 is
## NA, then size1/height1 is the only peak, so take that.

## Limit this to cases with at least 'min.comb' enzyme/primer
## combinations.

## Return a data.frame with sample.fk, primer, enzyme, size as
## columns.
find.max.peaks <- function(samples.data) {
  f <- function(idx) {
    i <- idx[order(samples.data$height[idx], decreasing=TRUE)[1:2]]
    c(size=samples.data$size[i], height=samples.data$height[i])
  }
  i <- paste(samples.data$sample.fk, samples.data$primer,
             samples.data$enzyme, sep="\r")
  out <- samples.data[tapply(1:nrow(samples.data), i, "[", 1),
                      c("sample.fk", "primer", "enzyme")]
  heights <- do.call(rbind, tapply(1:nrow(samples.data), i, f))
  colnames(heights)[c(1,3)] <- c("size", "height")
  out <- cbind(out, heights)
  out$ratio <- out$height/out$height2
  out[do.call(order, out),]
}


## (1) build.knowns: Automatically Build Knowns Database

## This function uses several fiters to select likely knowns.  Samples
## are considered useful if they have an adequate number of
## enzyme/primer combinations, and if for each combination they have
## either a single peak, or a peak that is *distinct enough* from all
## the others.

## For all samples and enzyme/primer combinations, the ratio of the
## largest to the second largest peak is calculated.  If it is greater
## than min.ratio, then that combination is accepted.  If the sample
## has at least min.comb valid enzyme/primer combinations, then that
## sample is included in the knowns database.  If min.comb is NA (the
## default), then every enzyme/primer combination present in the data
## is required.

## Most work is done by find.max.peaks()
## 'restrict' only allows use of cases where d$species is non-NA.
build.knowns <- function(d, min.ratio=3, min.comb=NA,
                         restrict=FALSE, ...) {
  valid.TRAMPsamples(d)
  samples.data <- d$data
  samples.info <- d$info

  if ( restrict ) {
    samples.info <- samples.info[!is.na(samples.info$species),]
    samples.data <- samples.data[samples.data$sample.fk %in%
                                 samples.info$sample.pk,]
  }
  
  if ( nrow(samples.data) < 1 )
    stop("No data to build knowns database from")

  if ( is.na(min.comb) )
    min.comb <- nrow(unique(samples.data[c("primer", "enzyme")]))

  knowns <- find.max.peaks(samples.data)
  knowns <- knowns[is.na(knowns$height2) | knowns$ratio > min.ratio,]

  knowns <- knowns[knowns$sample.fk %in%
                   names(which(table(knowns$sample.fk) >= min.comb)),]
  knowns <- knowns[c("sample.fk", "primer", "enzyme", "size")]
  names(knowns)[1] <- "knowns.fk"

  knowns.pk <- unique(knowns$knowns.fk)
  species <-
    samples.info$species[match(knowns.pk, samples.info$sample.pk)]
  knowns.info <- data.frame(knowns.pk=knowns.pk, species=species,
                            stringsAsFactors=FALSE)

  if ( nrow(knowns.info) == 0 )
    stop("Did not find any potential knowns (try a lower min.comb?)")
  
  TRAMPknowns(knowns, knowns.info, ...)
}

## _Basic_ indexing method.
"[.TRAMPknowns" <- function(x, i, na.interp=TRUE, ...) {
  valid.TRAMPknowns(x)
  if ( is.numeric(i) ) {
    if ( !all(i %in% labels(x)) )
      stop("Unknown knowns: ",
           paste(i[!(i %in% labels(x))], collapse=", "))
  } else if ( is.logical(i) ) {
    if ( length(i) != nrow(x$info) )
      stop("Logical index of incorrect length")
    i[is.na(i)] <- na.interp
    i <- x$info$knowns.pk[i]
  } else stop("Invalid index type")

  ## These generate NOTEs in R CMD CHECK
  x$info <- subset(x$info, x$info$knowns.pk %in% i)
  x$data <- subset(x$data, x$data$knowns.fk %in% i)
  ##x$info <- x$info[x$info$knowns.pk %in% i,]
  ##x$data <- x$data[x$data$knowns.pk %in% i,]
  tidy.TRAMPknowns(x)
}
