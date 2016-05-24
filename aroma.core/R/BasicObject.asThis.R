setMethodS3("asThis", "BasicObject", function(this, object, ...) {
  # Nothing to do?
  if (inherits(object, class(this)[1])) {
    return(object);
  }

  # Locate all asNnn() methods
  names <- sprintf("as%s", class(this));
  for (kk in seq_along(names)) {
    name <- names[kk];
    if (exists(name, mode="function")) {
      asFcn <- get(name, mode="function");
      res <- asFcn(object);
      return(res);
    }
  } # for (kk ...)

  throw("Could not coerce ", class(object)[1], " object to class ", 
        class(this)[1], ". Tried the following asNnn() methods: ", 
                                       paste(names, collapse=", "));
})

setMethodS3("asThis", "Object", function(this, object, ...) {
  asThis.BasicObject(this, object, ...);
})


############################################################################
# HISTORY:
# 2009-03-31
# o Created.
############################################################################
