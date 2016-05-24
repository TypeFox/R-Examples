setMethodS3("remap", "default", function(x, map, values=NULL, ...) {
  # Argument 'map':
  mode <- mode(x);
  if (!identical(mode(map), mode)) {
    throw("Argument 'map' is of a different mode than 'x': ",
                                        mode(map), " != ", mode);
  }

  # Argument 'values':
  if (is.null(values)) {
    values <- seq_along(map);
    mode(values) <- mode;
  }

  # Nothing todo?
  if (identical(map, values)) {
    return(x);
  }

  # Allocate return object
  y <- vector(mode(values), length=length(x));
  dim(y) <- dim(x);

  # Remap
  nbrOfValues <- length(map);
  for (kk in seq_len(nbrOfValues)) {
    idxs <- (x == map[kk]);
    idxs <- which(idxs);
    if (length(idxs) > 0L) {
      y[idxs] <- values[kk];
    }
  }

  y;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2008-07-09
# o Created.
############################################################################
