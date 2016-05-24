setMethodS3("getStateColorMap", "SegmentedGenomicSignalsInterface", function(this, ...) {
  colorMap <- getBasicField(this, ".stateColorMap");
  if (is.null(colorMap)) {
    this <- setStateColorMap(this, colorMap="default");
    colorMap <- getBasicField(this, ".stateColorMap");
  }
  colorMap;
})


setMethodS3("setStateColorMap", "SegmentedGenomicSignalsInterface", function(this, colorMap="default", ...) {
  # Argument 'colorMap':
  names <- names(colorMap);
  if (is.null(names)) {
    colorMap <- match.arg(colorMap);
    if (colorMap == "default") {
      colorMap <- c(
        "NA" = "#999999",  # missing values
        "0" = "#000000",   # neutral
        "-" = "blue",      # losses
        "+" = "red",       # gains
        "*"  = "#000000"   # default
      );
    }
  } else {
    colorMap <- Arguments$getCharacters(colorMap);
    names(colorMap) <- names;
  }

  this <- setBasicField(this, ".stateColorMap", colorMap);

  invisible(this);
})


setMethodS3("getStateColors", "SegmentedGenomicSignalsInterface", function(this, na.rm=FALSE, ...) {
  colorMap <- getStateColorMap(this);
  if (na.rm) {
    colorMap["NA"] <- as.character(NA);
  }
  hasDefColor <- is.element("*", names(colorMap));
  if (hasDefColor) {
    for (type in c("0", "-", "+")) {
      if (!is.element(type, names(colorMap))) {
        colorMap[type] <- colorMap["*"];
      }
    }
  }

  states <- getStates(this);
#  print(table(states, exclude=NULL));
  uStates <- sort(unique(states), na.last=TRUE);
  uStates <- na.omit(uStates);

  # Default missing-value colors
  naColor <- as.character(colorMap["NA"]);
  cols <- rep(naColor, times=length(states));

  # Neutral states
  if (is.element("0", names(colorMap))) {
    idxs <- which(states == 0);
    cols[idxs] <- colorMap["0"];
  }

  # Negative states
  if (is.element("-", names(colorMap))) {
    idxs <- which(states < 0);
    cols[idxs] <- colorMap["-"];
  }

  # Positive states
  if (is.element("+", names(colorMap))) {
    idxs <- which(states > 0);
    cols[idxs] <- colorMap["+"];
  }

  for (kk in seq_along(uStates)) {
    state <- uStates[kk];
    key <- sprintf("%s", state);

    if (is.element(key, names(colorMap))) {
      idxs <- which(states == state);
      cols[idxs] <- colorMap[key];
    }
  } # for (kk ...)

  cols;
})


############################################################################
# HISTORY:
# 2009-07-01
# o BUG FIX: getStateColorMap() with "*" gave all NA colors.
# 2009-06-29
# o Created.
############################################################################
