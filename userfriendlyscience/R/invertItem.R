### To invert mirrored items
invertItem <- function(item, range=NULL, ignorePreviousInversion = FALSE) {
  ### Check whether this was already inverted
  if (!is.null(attr(item, "inverted"))) {
    if ((attr(item, "inverted") == TRUE) & !(ignorePreviousInversion)) {
      warning("Vector '", substitute(deparse(item)),
              "' has already been inverted! ",
              "Set ignorePreviousInversion to TRUE to override this ",
              "check and invert the vector anyway.");
    }
  }
  
  ### Not inverted yet (or ignorePreviousInversion set to TRUE)
  if (is.numeric(item)) {
    if (is.null(range)) {
      res <- sum(range(na.omit(item))) - item;
    }
    else {
      res <- sum(range(range)) - item;
    }
  }
  else {
    stop("Provide a numeric vector!");
  }
  attr(res, "inverted") <- TRUE;
  return(res);
}