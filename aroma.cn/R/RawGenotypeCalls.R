setConstructorS3("RawGenotypeCalls", function(...) {
  extend(RawAlleleBFractions(...), "RawGenotypeCalls")
})

setMethodS3("getCalls", "RawGenotypeCalls", function(this, flavor=c("fracB", "AB"), ...) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  y <- getSignals(this);
  if (flavor == "fracB") {
    calls <- y;
  } else if (flavor == "AB") {
    naValue <- as.character(NA);
    calls <- rep(naValue, times=nbrOfLoci(this));
    calls[y == 0] <- "AA";
    calls[y == 1/2] <- "AB";
    calls[y == 1] <- "BB";
  }
  calls;
})

setMethodS3("isHeterozygous", "RawGenotypeCalls", function(this, ...) {
  calls <- getCalls(this, flavor="fracB", ...);
  res <- (is.finite(calls) & (calls == 1/2));
  res;
})

setMethodS3("isHomozygous", "RawGenotypeCalls", function(this, ...) {
  calls <- getCalls(this, flavor="fracB", ...);
  res <- (is.finite(calls) & (calls != 1/2));
  res;
})

setMethodS3("getColors", "RawGenotypeCalls", function(this, colorMap=c(AA="red", AB="black", BB="red", "NA"="#999999"), ...) {
  # Argument 'colorMap':
  colorMap <- Arguments$getCharacters(colorMap, useNames=TRUE);
  calls <- getCalls(this, flavor="AB");
  col <- colorMap[calls];
  col[is.na(col)] <- colorMap["NA"];
  col;
})



############################################################################
# HISTORY:
# 2009-10-10
# o Added RawGenotypeCalls.
# o Created.
############################################################################ 
