setMethodS3("getValueAt", "CopyNumberRegions", function(this, field, at, ...) {
  df <- as.data.frame(this);
  values <- df[[field]];
  idx <- which(df$start <= at & at < df$stop);
  if (length(idx) == 1) {
    res <- values[idx];
  } else {
    res <- as.double(NA);
  }
  res;
}, protected=TRUE)


setMethodS3("getMeanAt", "CopyNumberRegions", function(this, ...) {
  getValueAt(this, field="mean", ...);
})


setMethodS3("prune", "CopyNumberRegions", function(this, delta=0, ...) {
  # Argument 'delta':
  delta <- Arguments$getDouble(delta);

  data <- as.data.frame(this);
  rr <- 1L;
  while (rr < nrow(data)) {
    regionA <- data[rr,];
    regionB <- data[rr+1,];
    isEqual <- (abs(regionA$mean - regionB$mean) <= delta)
    if (isEqual) {
      data[rr,"start"] <- min(regionA$start, regionB$start);
      data[rr,"stop"] <- max(regionA$stop, regionB$stop);
      data[rr,"count"] <- sum(c(regionA$count, regionB$count), na.rm=TRUE);
      data <- data[-(rr+1),,]
    } else {
      rr <- rr + 1L;
    }
  }

  res <- clone(this);
  for (field in colnames(data)) {
    res[[field]] <- data[[field]];
  }

  res;
})


setMethodS3("+", "CopyNumberRegions", function(e1, e2) {
  # To please R CMD check
  this <- e1;
  other <- e2;

  # Argument 'other':
  other <- Arguments$getInstanceOf(other, class(this)[1]);

  chromosome <- this$chromosome;
  if (length(unique(chromosome)) > 1) {
    throw("Adding regions across multiple chromosomes is not supported.");
  }
  chromosome <- chromosome[1];

  if (!identical(unique(other$chromosome), chromosome)) {
    throw("Argument 'other' is for a different chromosome: ", 
          unique(other$chromosome), " != ", chromosome);
  }

  regList <- list(this, other);
  dfList <- lapply(regList, FUN=as.data.frame);

  # Identify all unique change points
  cps <- lapply(dfList, FUN=function(df) c(df[,c("start","stop")]));
  cps <- unlist(cps, use.names=FALSE);
  cps <- unique(sort(cps));

  starts <- cps[-length(cps)];
  stops <- cps[-1];

  res <- clone(this);
  res$chromosome <- rep(chromosome, times=length(starts));
  res$start <- starts;
  res$stop <- stops;

  fields <- c("mean");
  for (field in fields) {
    values <- double(length(starts));

    for (rr in seq_along(starts)) {
      region <- c(starts[rr], stops[rr]);
      mu <- sapply(regList, FUN=function(reg) {
        getValueAt(reg, field=field, at=region[1]);
      });
      mu <- unlist(mu, use.names=FALSE);
      mu <- sum(mu, na.rm=TRUE);
      values[rr] <- mu;
    }
    res[[field]] <- values;
  }

##  # Merge equal regions
##  res <- prune(res, delta=0);

  res;  
}, appendVarArgs=FALSE, validators=NULL)


setMethodS3("-", "CopyNumberRegions", function(e1, e2) {
  # To please R CMD check
  this <- e1;
  other <- e2;

  otherNeg <- clone(other);
  otherNeg$mean <- -otherNeg$mean;
  res <- this + otherNeg;
  res;
}, appendVarArgs=FALSE, validators=NULL)


setMethodS3("*", "CopyNumberRegions", function(e1, e2) {
  # To please R CMD check
  this <- e1;
  value <- e2;

  # Swap 'this' and 'value'?
  if (inherits(value, "CopyNumberRegions")) {
    tmp <- this;
    this <- value;
    value <- tmp;
  }

  value <- Arguments$getDouble(value);

  res <- clone(this);
  res$mean <- value*res$mean;

  res;
}, appendVarArgs=FALSE, validators=NULL)


setMethodS3("xRange", "CopyNumberRegions", function(this, ...) {
  range(c(this$start, this$stop), na.rm=TRUE);
})

setMethodS3("xMin", "CopyNumberRegions", function(this, ...) {
  xRange(this)[1];
})

setMethodS3("xMax", "CopyNumberRegions", function(this, ...) {
  xRange(this)[2];
})


setMethodS3("simulateRawCopyNumbers", "CopyNumberRegions", function(this, x=seq(from=xMin(this), to=xMax(this), length=n), n, rfcn=function(n, x, mu, ...) rnorm(n=n, mean=0, sd=sd), sd=0.1, ...) {
  # Argument 'x':
  x <- Arguments$getNumerics(x);

  # Argument 'rfcn':
  if (!is.function(rfcn)) {
    throw("Argument 'rfcn' is not a function: ", class(rfcn)[1]);
  }


  # Allocate result
  nbrOfLoci <- length(x);
  naValue <- as.double(NA);
  mu <- rep(naValue, times=nbrOfLoci);

  data <- as.data.frame(this);
  for (rr in seq_len(nrow(data))) {
    region <- c(data[rr,"start"], data[rr,"stop"]);
    muRR <- data[rr,"mean"];
    idxs <- which(region[1] <= x & x < region[2]);
    if (length(idxs) > 0) {
      mu[idxs] <- muRR;
    }
  } # for (rr ...)

  eps <- rfcn(n=n, x=x, mu=mu);
  # Sanity check
  stopifnot(length(eps) == nbrOfLoci);

  y <- mu + eps;

  # Sanity check
  stopifnot(length(y) == nbrOfLoci);

  chromosome <- this$chromosome[1];
  res <- RawCopyNumbers(x=x, cn=y, chromosome=chromosome);

  res;
}, protected=TRUE)


## setMethodS3("simulateGaussian", "CopyNumberRegions", function(this, ..., sd=1) {
##  simulate(this, ..., rfcn=function(n, ...) rnorm(n=n, mean=0, sd=sd));
## })


############################################################################
# HISTORY:
# 2011-12-10
# o ROBUSTNESS: Turned of RCC validation for "+", "-" and "*" methods.
# 2010-07-24
# o Added simulateRawCopyNumbers() for CopyNumberRegions.
# o Added xRange(), xMin() and xMax() for CopyNumberRegions.
# o Added "*"() to multiple a CopyNumberRegions profile with a scalar.
# o Added "+"() and "-"() to add two CopyNumberRegions profiles.
# o Added getMeanAt().
# o Added prune().
############################################################################
