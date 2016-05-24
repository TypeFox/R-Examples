setConstructorS3("CopyNumberRegions", function(chromosome=NULL, start=NULL, stop=NULL, mean=NULL, count=NULL, call=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(chromosome)) {
    # Argument 'start':
    start <- Arguments$getNumerics(start);
    n <- length(start);

    # Arguments 'stop' & 'mean':
    stop <- Arguments$getNumerics(stop, length=c(n,n));
    mean <- Arguments$getNumerics(mean, length=c(n,n));

    # Argument 'chromosome':
    if (length(chromosome) == 1) {
      chromosome <- rep(chromosome, times=n);
    }
    chromosome <- Arguments$getIntegers(chromosome, length=c(n,n));

    # Argument 'count':
    if (!is.null(count)) {
      count <- Arguments$getIntegers(count, length=c(n,n));
    }
  }

  extend(Object(), "CopyNumberRegions",
    chromosome = chromosome,
    start = start,
    stop = stop,
    mean = mean,
    count = count,
    call= call,
    ...
  )
})


setMethodS3("as.character", "CopyNumberRegions", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Number of regions: %d", nbrOfRegions(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)


setMethodS3("nbrOfRegions", "CopyNumberRegions", function(this, ...) {
  length(this$start);
})

setMethodS3("getLength", "CopyNumberRegions", function(this, ...) {
  this$stop - this$start;
})


setMethodS3("getDensity", "CopyNumberRegions", function(this, field="mean", adjust=0.2, ...) {
  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0,Inf));

  # Region signal
  y <- this[[field]];

  # Weight by length of regions
  w <- getLength(this);

  # Drop missing values
  keep <- is.finite(y);
  y <- y[keep];
  w <- w[keep];

  # Special case
  if (length(y) == 1) {
    y <- rep(y, times=2);
    w <- rep(w, times=2);
  }

  # Standardize weights
  w <- w / sum(w, na.rm=TRUE);

  # Estimate the weighted empircal density
  density(y, weights=w, adjust=adjust, ...);
})


setMethodS3("as.data.frame", "CopyNumberRegions", function(x, ...) {
  # To please R CMD check
  this <- x;

  fields <- c("chromosome", "start", "stop", "mean", "count", "call");
  data <- lapply(fields, FUN=function(field) this[[field]]);
  names(data) <- fields;
  data <- data[!sapply(data, is.null)];
  data <- as.data.frame(data);
  data;
})


setMethodS3("extractIGV", "CopyNumberRegions", function(this, ...) {
  data <- as.data.frame(this, ...);
  names <- colnames(data);
  names <- gsub("chromosome", "Chromosome", names);
  names <- gsub("start", "Start Position", names);
  names <- gsub("stop", "End Position", names);
  names <- gsub("count", "Num markers", names);
  names <- gsub("mean", "Seg CN", names);
  colnames(data) <- names;

  # Append 'Sample' column
  name <- this$name;
  if (is.null(name)) {
    name <- "<Unknown Sample>";
  }
  data <- cbind(Sample=name, data);

  # Reorder for IGV
  cols <- c("Sample", "Chromosome", "Start Position",
            "End Position", "Num markers", "Seg CN");
  data <- data[,cols,drop=FALSE];

  data;
})


setMethodS3("equals", "CopyNumberRegions", function(this, other, ...) {
  dfThis <- as.data.frame(this);
  dfOther <- as.data.frame(other);
  res <- all.equal(dfThis, dfOther);
  isTRUE(res);
})


setMethodS3("applyRows", "CopyNumberRegions", function(this, FUN, ...) {
  data <- as.data.frame(this);
  res <- vector("list", nrow(data));

  if (nrow(data) > 0) {
    o <- order(data[,"chromosome"], data[,"start"]);
    data <- data[o,,drop=FALSE];
    for (kk in seq_len(nrow(data))) {
      res[[kk]] <- FUN(data[kk,,drop=FALSE], ...);
    }
  }

  res;
}, protected=TRUE)


setMethodS3("drawLevels", "CopyNumberRegions", function(this, col="red", lwd=2, lty=1, xScale=1e-6, yScale=1, ...) {
  col0 <- col;
  lwd0 <- lwd;
  lty0 <- lty;
  res <- applyRows(this, FUN=function(cnr) {
    x <- c(cnr[["start"]], cnr[["stop"]]);
    y <- rep(cnr[["mean"]], times=2);
    if (is.function(col0))
      col <- col0(cnr);
    if (is.function(lwd0))
      lwd <- lwd0(cnr);
    if (is.function(lty0))
      lty <- lty0(cnr);
    lines(x=xScale*x, y=yScale*y, col=col, lwd=lwd, lty=lty, ...);
  });
  invisible(res);
})



setMethodS3("lines", "CopyNumberRegions", function(x, col="red", lwd=2, xScale=1e-6, yScale=1, ...) {
  # To please R CMD check.
  this <- x;

  data <- as.data.frame(this);
  if (nrow(data) > 0) {
    o <- order(data[,"start"]);
    data <- data[o,,drop=FALSE];
    xx <- t(data[,c("start", "stop"),drop=FALSE]);
    yy <- rep(this$mean[o], each=2);
    lines(x=xScale*xx, y=yScale*yy, col=col, lwd=lwd, ...);
  }
})



setMethodS3("subset", "CopyNumberRegions", function(x, subset, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'subset':
  if (all(subset < 0)) {
    subset <- -Arguments$getIndices(-subset, max=nbrOfRegions(this));
  } else {
    subset <- Arguments$getIndices(subset, max=nbrOfRegions(this));
  }

  # Get field names.  AD HOC /HB 2010-07-19
  data <- as.data.frame(this);
  fields <- colnames(data);
  # Not needed anymore
  data <- NULL;

  res <- clone(this);
  for (field in fields) {
    res[[field]] <- res[[field]][subset];
  }

  res;
})


setMethodS3("extractCopyNumberRegions", "default", abstract=TRUE);

setMethodS3("extractCNRs", "default", function(...) {
  extractCopyNumberRegions(...);
})


############################################################################
# HISTORY:
# 2010-09-12
# o Added getDensity() for CopyNumberRegions.
# o Added virtual field 'length'/getLength() for CopyNumberRegions.
# 2010-07-19
# o Added subset() for CopyNumberRegions.
# 2010-04-06
# o Added equals() for CopyNumberRegions.
# 2009-05-16
# o Now the constructor CopyNumberRegions() coerce numerics only if
#   necessary, i.e. it keeps integers if integers, otherwise to doubles.
#   This is a general design of aroma.* that saves some memory.
# 2009-05-13
# o CLEAN UP: drawLevels() no longer returns a (useless) list structure.
# 2008-05-17
# o Added abstract default extractCopyNumberRegions().
# o Moved to aroma.core.
# 2008-04-17
# o BUG FIX: applyRows() and lines() of CopyNumberRegions did not handle
#   cases with zero regions.
# 2007-09-04
# o Now CopyNumberRegions also contains an 'chromosome' field.
# o BUG FIX: as.data.frame() gave an error if some optional fields were
#   NULL.
# 2007-08-22
# o Created.  Need a generic container for holding copy number regions and
#   to plot them nicely.
############################################################################
