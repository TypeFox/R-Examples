###########################################################################/**
# @RdocClass AbstractCBS
#
# @title "The AbstractCBS class"
#
# \description{
#  @classhierarchy
#
#  All CBS-style segmentation results extend this class, e.g.
#  @see "CBS" and @see "PairedPSCBS".
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A @list structure containing the segmentation results.}
#   \item{sampleName}{A @character string.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setConstructorS3("AbstractCBS", function(fit=list(), sampleName=fit$sampleName, ...) {
  # Argument 'sampleName':
  if (!is.null(sampleName)) {
    sampleName <- Arguments$getCharacter(sampleName);
  }

  fit$sampleName <- sampleName;
  extend(fit, "AbstractCBS");
})


setMethodS3("print", "AbstractCBS", function(x, ...) {
  # To please R CMD check
  fit <- x;

  segs <- getSegments(fit, simplify=TRUE, ...);
  print(segs);
}, protected=TRUE)




setMethodS3("all.equal", "AbstractCBS", function(target, current, check.attributes=FALSE, ...) {
  # NOTE: Here we cannot trust argument '...', because it may contain
  # copies of 'target' and 'current'
  args <- list(...);
  drop <- integer(0L);
  for (kk in seq_along(args)) {
    if (identical(args[[kk]], target)) drop <- c(drop, kk);
    if (identical(args[[kk]], current)) drop <- c(drop, kk);
  }
  if (length(drop) > 0L) {
    args <- args[-drop];
    str(args);
#    assign("...", args, inherits=FALSE);
  }
  args <- list(...);

  # Compare class attributes
  res <- all.equal(class(target), class(current));
  if (!isTRUE(res)) {
    return(res);
  }

  # Compare locus-level data
  dataT <- getLocusData(target);
  dataC <- getLocusData(current);
  res <- all.equal(dataT, dataC, check.attributes=check.attributes);
  if (!isTRUE(res)) {
    attr(res, "what") <- "getLocusData()";
    return(res);
  }

  # Compare segments
  dataT <- getSegments(target);
  dataC <- getSegments(current);
  res <- all.equal(dataT, dataC, check.attributes=check.attributes);
  if (!isTRUE(res)) {
    attr(res, "what") <- "getSegments()";
    return(res);
  }

  # Compare field names
  fieldsT <- names(target);
  fieldsC <- names(current);
  res <- all.equal(fieldsT, fieldsC, check.attributes=check.attributes);
  if (!isTRUE(res)) {
      attr(res, "what") <- "names";
    return(res);
  }

  # Compare other fields
  for (key in fieldsT) {
    dataT <- target[[key]];
    dataC <- current[[key]];
    res <- all.equal(dataT, dataC, check.attributes=check.attributes);
    if (!isTRUE(res)) {
      attr(res, "what") <- sprintf("[[\"%s\"]]", key);
      return(res);
    }
  } # for (key ...)

  return(TRUE);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod save
#
# @title "Saves an AbstractCBS object to file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to @see "R.utils::saveObject".}
# }
#
# \value{
#   Returns what @see "R.utils::saveObject" returns.
# }
#
# @author
#
# \seealso{
#   Internally @see "R.utils::saveObject" is used.
#   To load an object, see @seemethod "load".
#   @seeclass.
# }
#*/###########################################################################
setMethodS3("save", "AbstractCBS", function(this, ...) {
  saveObject(this, ...);
})


###########################################################################/**
# @RdocMethod load
# @alias load
#
# @title "Loads an AbstractCBS object from file"
#
# \description{
#  @get "title" and assert that it is of the requested class.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to @see "R.utils::loadObject".}
# }
#
# \value{
#   Returns the loaded AbstractCBS object.
# }
#
# @author
#
# \seealso{
#   Internally @see "R.utils::loadObject" is used.
#   To save an object, see @seemethod "save".
#   @seeclass.
# }
#*/###########################################################################
setMethodS3("load", "AbstractCBS", function(static, ...) {
  object <- loadObject(...);

  # Patch for changes in class structure in PSCBS v0.13.2 -> v0.13.3.
  if (!inherits(object, "AbstractCBS")) {
    if (inherits(object, "CBS")) {
      class(object) <- c(class(object), "AbstractCBS");
      warning("Added 'AbstractCBS' to the class hierarchy of the loaded ", class(object)[1], " object.");
    } else if (inherits(object, "PairedPSCBS")) {
      class(object) <- c(class(object), "AbstractCBS");
      warning("Added 'AbstractCBS' to the class hierarchy of the loaded ", class(object)[1], " object.");
    }
  }

  # Sanity check
  if (!inherits(object, class(static)[1])) {
    throw("Loaded an object from file, but it does not inherit from ",
          class(static)[1], " as expected: ", hpaste(class(object)));
  }

  object;
}, static=TRUE)



###########################################################################/**
# @RdocMethod getSampleName
# @aliasmethod sampleName
#
# @title "Gets the name of the sample segmented"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "setSampleName".
#   @seeclass.
# }
#*/###########################################################################
setMethodS3("getSampleName", "AbstractCBS", function(fit, ...) {
  name <- fit$sampleName;
  if (is.null(name)) {
    name <- as.character(NA);
  }
  name;
}, protected=TRUE)

setMethodS3("sampleName", "AbstractCBS", function(fit, ...) {
  getSampleName(fit);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod setSampleName
# @aliasmethod sampleName<-
#
# @title "Sets the name of the sample segmented"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{A @character string.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) an updated object.
# }
#
# @author
#
# \seealso{
#   @seeclass.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("setSampleName", "AbstractCBS", function(fit, name, ...) {
  # Argument 'value':
  name <- Arguments$getCharacter(name);

  fit$sampleName <- name;

  invisible(fit);
}, protected=TRUE)


setMethodS3("sampleName<-", "AbstractCBS", function(x, value) {
  setSampleName(x, value);
}, protected=TRUE, addVarArgs=FALSE)

"sampleName<-" <- function(x, value) {
  UseMethod("sampleName<-");
}



###########################################################################/**
# @RdocMethod getLocusData
# @aliasmethod setLocusData
# @alias setLocusData.AbstractCBS
#
# @title "Gets the locus-level data"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{splitters}{If @TRUE, "splitters" between chromosomes are
#     preserved, otherwise dropped.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a JxL @data.frame, where J in the number of loci,
#   and L is the number of locus-specific fields.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getLocusData", "AbstractCBS", abstract=TRUE);

setMethodS3("setLocusData", "AbstractCBS", function(fit, loci, ...) {
  # Argument 'loci':
  loci <- Arguments$getInstanceOf(loci, "data.frame");
  nbrOfLoci <- nbrOfLoci(fit);
  if (nrow(loci) != nbrOfLoci) {
    throw("Cannot set locus-level data. The number of loci to be set differ from the existing number of loci: ", nrow(loci), " != ", nbrOfLoci);
  }

  fit$data <- loci;

  invisible(fit);
}, protected=TRUE)


setMethodS3("getLocusSignalNames", "AbstractCBS", abstract=TRUE, protected=TRUE);

setMethodS3("getSegmentTrackPrefixes", "AbstractCBS", abstract=TRUE, protected=TRUE);


###########################################################################/**
# @RdocMethod nbrOfLoci
#
# @title "Gets the number of loci"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{splitters, ...}{Arguments passed to @seemethod "getLocusData".}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfLoci", "AbstractCBS", function(fit, splitters=FALSE, ...) {
  data <- getLocusData(fit, splitters=splitters, ...);
  nrow(data);
})



###########################################################################/**
# @RdocMethod getSegments
# @aliasmethod setSegments
# @alias setSegments.AbstractCBS
#
# @title "Gets the segments"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{simplify}{If @TRUE, redundant and intermediate information is dropped.}
#  \item{splitters}{If @TRUE, "splitters" between chromosomes are
#     preserved, otherwise dropped.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a SxK @data.frame, where S in the number of segments,
#   and K is the number of segment-specific fields.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getSegments", "AbstractCBS", abstract=TRUE);


setMethodS3("setSegments", "AbstractCBS", function(fit, segments, ...) {
  # Argument 'segments':
  segments <- Arguments$getInstanceOf(segments, "data.frame")
  nbrOfSegs <- nbrOfSegments(fit, ...)
  if (nrow(segments) != nbrOfSegs) {
    throw("Cannot set segments. The number of segments to be set differ from the existing number of segments: ", nrow(segments), " != ", nbrOfSegs)
  }

  fit$output <- segments

  invisible(fit)
}, protected=TRUE)


setMethodS3("getChangePoints", "AbstractCBS", abstract=TRUE);



###########################################################################/**
# @RdocMethod resetSegments
#
# @title "Reset the segments"
#
# \description{
#   @get "title".  More precisely, it removes columns in the segmentation
#   result table that have been added by methods after the actual
#   segmentation method, e.g. bootstrap estimated mean level quantiles
#   and various calls.
#   It leave the basic segmentation results untouched,
#   i.e. the partitioning and the segment means.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an object if the same class as the input result.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("resetSegments", "AbstractCBS", function(fit, ...) {
  segs <- getSegments(fit, splitters=TRUE)
  names <- colnames(segs);

  excl <- NULL;

  # Drop all quantile mean level estimates (from bootstrapping)
  idxs <- grep("_[0-9.]*[%]$", names);
  excl <- c(excl, idxs);

  # Drop all calls
  idxs <- grep("Call$", names);
  excl <- c(excl, idxs);

  excl <- unique(excl);
  if (length(excl) > 0L) {
    segs <- segs[,-excl];
  }

  fit <- setSegments(fit, segs, splitters=TRUE)
  invisible(fit);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod nbrOfSegments
#
# @title "Gets the number of segments"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{splitters, ...}{Arguments passed to @seemethod "getSegments".}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfChangePoints"
#   @seemethod "nbrOfChromosomes"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfSegments", "AbstractCBS", function(this, splitters=FALSE, ...) {
  nrow(getSegments(this, splitters=splitters, ...));
})



###########################################################################/**
# @RdocMethod nbrOfChangePoints
#
# @title "Gets the number of change points"
#
# \description{
#   @get "title", which is defined as the number of segments minus
#   the number of chromosomes.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfSegments"
#   @seemethod "nbrOfChromosomes"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfChangePoints", "AbstractCBS", function(fit, ignoreGaps=FALSE, dropEmptySegments=TRUE, ...) {
  segs <- getSegments(fit, splitters=TRUE, addGaps=!ignoreGaps);
  if (dropEmptySegments) {
    prefix <- getSegmentTrackPrefixes(fit);
    keys <- sapply(prefix, FUN=function(x) {
      toCamelCase(paste(c(x, "nbr of loci"), collapse=" "));
    });
    counts <- as.matrix(segs[,keys]);
    counts <- rowSums(counts, na.rm=TRUE);
    segs$chromosome[counts == 0L] <- NA;
  }
  sum(!is.na(diff(segs$chromosome)));
})



###########################################################################/**
# @RdocMethod as.data.frame
#
# @title "Gets the table of segments"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame, where each row corresponds to
#   a unique segment.
# }
#
# @author
#
# \seealso{
#   Utilizes @seemethod "getSegments".
#   @seeclass.
# }
#*/###########################################################################
setMethodS3("as.data.frame", "AbstractCBS", function(x, ...) {
  getSegments(x, ...);
}, protected=TRUE)




###########################################################################/**
# @RdocMethod getChromosomes
#
# @title "Gets the set of chromosomes"
#
# \description{
#   @get "title" in the segmentation result.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getSegments".}
# }
#
# \value{
#   Returns a unique and sorted @vector of chromosomes segmented.
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfChromosomes".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChromosomes", "AbstractCBS", function(this, ...) {
  segs <- getSegments(this, ...);
  chromosomes <- sort(unique(segs$chromosome), na.last=TRUE);

  # Drop NA dividers
  if (length(chromosomes) > 1) {
    chromosomes <- chromosomes[!is.na(chromosomes)];
  }

  chromosomes;
})


###########################################################################/**
# @RdocMethod nbrOfChromosomes
#
# @title "Gets the number of chromosomes"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getChromosomes".}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seemethod "getChromosomes".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfChromosomes", "AbstractCBS", function(this, ...) {
  length(getChromosomes(this, ...));
})


setMethodS3("getSegmentSizes", "AbstractCBS", function(fit, by=c("length", "count"), ...) {
  by <- match.arg(by);

  if (by == "length") {
    prefix <- getSegmentTrackPrefixes(fit)[1];
    keys <- toCamelCase(paste(prefix, " ", c("start", "end")));
  } else if (by == "count") {
    keys <- "nbrOfLoci";
  }
  data <- getSegments(fit, ...)[,keys];

  if (by == "length") {
    res <- data[[2L]]-data[[1L]]+1L;
  } else if (by == "count") {
    res <- data[[1L]];
  }
  res;
})


setMethodS3("extractCNs", "AbstractCBS", abstract=TRUE);

setMethodS3("sampleCNs", "AbstractCBS", function(fit, size=NULL, ...) {
  data <- extractCNs(fit, ...);

  if (!is.null(size)) {
    sizes <- getSegmentSizes(fit, ...);
    # Sanity check
    stopifnot(length(sizes) == nrow(data));
    idxs <- sample(nrow(data), size=size, replace=TRUE, prob=sizes);
    data <- data[idxs,,drop=FALSE];
  }

  data;
})

###########################################################################/**
# @RdocMethod updateMeans
# @alias updateMeans.CBS
# @alias updateMeans.NonPairedPSCBS
# @alias updateMeans.PairedPSCBS
#
# @title "Updates the CN mean levels for each segment independently"
#
# \description{
#  @get "title" as if they were one large segment.
#  The locus-level data is not updated/modified.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments specific to the class.}
# }
#
# \value{
#   Returns an object of the same class.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("updateMeans", "AbstractCBS", abstract=TRUE, protected=TRUE);


setMethodS3("getMeanEstimators", "AbstractCBS", function(fit, which=NULL, default=mean, ...) {
  estList <- fit$params$meanEstimators;
  if (is.null(estList)) {
    estList <- list();
  }

  if (is.null(which)) which <- names(estList);

  for (key in which) {
    fcn <- estList[[key]];
    if (is.null(fcn)) {
      fcn <- default;
    } else if (is.character(fcn)) {
      fcn <- get(fcn, mode="function");
    }
    estList[[key]] <- fcn;
  }

  estList;
}, protected=TRUE)


setMethodS3("setMeanEstimators", "AbstractCBS", function(fit, ...) {
  estList <- fit$params$meanEstimators;
  if (is.null(estList)) {
    estList <- list();
  }

  args <- list(...);

  # Nothing todo?
  if (length(args) == 0L) {
    return(invisible(fit));
  }

  keys <- names(args);
  if (is.null(keys)) {
    throw("Estimators arguments must be named.");
  }

  for (key in keys) {
    fcn <- args[[key]];
    if (is.function(fcn)) {
    } else if (is.character(fcn)) {
      if (!exists(fcn, mode="function")) {
        throw(sprintf("No such '%s' estimator function: %s", key, fcn));
      }
    } else {
      throw(sprintf("Estimator argument '%s' must be a function or character string: %s", key, mode(fcn)));
    }
    estList[[key]] <- fcn;
  }

  fit$params$meanEstimators <- estList;

  invisible(fit);
}, protected=TRUE)


setMethodS3("resegment", "AbstractCBS", abstract=TRUE, protected=TRUE);


setMethodS3("getChromosomeRanges", "AbstractCBS", abstract=TRUE, protected=TRUE);

setMethodS3("getChromosomeOffsets", "AbstractCBS", function(fit, resolution=1e6, ...) {
  # Argument 'resolution':
  if (!is.null(resolution)) {
    resolution <- Arguments$getDouble(resolution, range=c(1,Inf));
  }

  data <- getChromosomeRanges(fit, ...);
  splits <- data[,"start"] + data[,"length"];

  if (!is.null(resolution)) {
    splits <- ceiling(splits / resolution);
    splits <- resolution * splits;
  }

  offsets <- c(0L, cumsum(splits));
  names(offsets) <- c(rownames(data), NA);

  offsets;
}, protected=TRUE) # getChromosomeOffsets()



###########################################################################/**
# @RdocMethod ploidy
# @aliasmethod ploidy<-
# @aliasmethod setPloidy
# @aliasmethod adjustPloidyScale
# @alias adjustPloidyScale.PairedPSCBS
# @alias adjustPloidyScale
# @alias ploidy
# @alias ploidy<-
# @alias setPloidy
#
# @title "Gets and sets ploidy"
#
# \description{
#  @get "title".
# }
#
# \usage{
#   \method{ploidy}{AbstractCBS}(fit, ...)
#   \method{ploidy}{AbstractCBS}(fit) <- value
# }
#
# \arguments{
#   \item{fit}{An @see "AbstractCBS" object.}
#   \item{value}{An @integer (in \eqn{1,2,\ldots}) specifying the genome ploidy .}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) an updated object.
# }
#
# @author
#
# \seealso{
#   @seeclass.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("ploidy", "AbstractCBS", function(fit, ...) {
  ploidy <- fit$params$ploidy;
  if (is.null(ploidy)) ploidy <- 2L;
  ploidy;
})

setMethodS3("ploidy<-", "AbstractCBS", function(fit, value) {
  fit <- setPloidy(fit, ploidy=value, update=TRUE);
  invisible(fit);
})

"ploidy<-" <- function(fit, value) {
  UseMethod("ploidy<-");
}

setMethodS3("setPloidy", "AbstractCBS", function(fit, ploidy=2L, update=TRUE, ...) {
  # Argument 'ploidy':
  ploidy <- Arguments$getInteger(ploidy, range=c(1,Inf));

  if (update) {
    # Calculate rescaling factor
    oldPloidy <- ploidy(fit);
    scale <- ploidy / oldPloidy;

    # Nothing todo?
    if (scale != 1) {
      fit <- adjustPloidyScale(fit, scale=scale, ...);
    }
  }

  fit$params$ploidy <- ploidy;
  invisible(fit);
}, protected=TRUE)


setMethodS3("adjustPloidyScale", "AbstractCBS", abstract=TRUE);


############################################################################
# HISTORY:
# 2013-11-05
# o Added basic implementations of setLocusData() and setSegments()
#   for AbstractCBS.
# 2013-10-20
# o Added abstract getChangePoints().
# 2013-05-07
# o Added ploidy() and ploidy()<- for AbstractCBS.
# 2013-02-01
# o Added resetSegments() for AbstractCBS, which drops extra segments
#   columns (e.g. bootstrap statisistics and calls) except those
#   obtained from the segment algorithm.
# 2013-01-15
# o Added get-/setMeanEstimators() for AbstractCBS.
# 2012-09-21
# o Now nbrOfChangePoints() for AbstractCBS calculates only change points
#   of connected neighboring segments.
# 2012-09-14
# o GENERALIZATION: Added getSegmentSizes() for AbstractCBS.
# o GENERALIZATION: Added getChromosomeOffsets() for AbstractCBS.
# 2012-09-13
# o Updated all.equal() for AbstractCBS to compare locus-level data,
#   segments, and other fields.
# 2012-06-03
# o DOCUMENTATION: Added Rd help for updateMeans().
# 2011-12-03
# o Now print() for AbstractCBS returns getSegments(..., simplify=TRUE).
# 2011-11-17
# o Added resegment() for AbstractCBS.
# 2011-10-30
# o Added save() and load() methods to AbstractCBS.
# 2011-10-16
# o Added sampleCNs() for AbstractCBS.
# o Added abstract getSegmentSizes() for AbstractCBS.
# o Added abstract extractCNs() for AbstractCBS.
# 2011-10-08
# o Added abstract updateMeans() for AbstractCBS.
# o Added all.equal() for AbstractCBS.
# o Added nbrOfChangePoints() for AbstractCBS.
# 2011-10-02
# o Created.
############################################################################
