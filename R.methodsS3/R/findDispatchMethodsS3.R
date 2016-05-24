###########################################################################/**
# @RdocDefault findDispatchMethodsS3
#
# @title "Finds the S3 methods that a generic function would call"
#
# \description{
#  @get "title", ordered according to an S3 @see "base::class" @vector.
# }
#
# @synopsis
#
# \arguments{
#   \item{methodName}{A @character string specifying the name of a
#     generic function.}
#   \item{classNames}{A @character @vector of @see "base::class" names.}
#   \item{firstOnly}{If @TRUE, only the first method is returned.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a names @list structure.
# }
#
# \seealso{
#   @see "getDispatchMethodS3".
# }
#
# @author
#
# @keyword programming
# @keyword methods
# @keyword internal
#*/###########################################################################
setMethodS3("findDispatchMethodsS3", "default", function(methodName, classNames, firstOnly=FALSE, ...) {
  # Argument 'methodName':
  methodName <- as.character(methodName);
  if (length(methodName) == 0) {
    throw("Argument 'methodName' is empty.");
  }
  if (length(methodName) > 1) {
    throw("Argument 'methodName' must only contain one element: ", paste(head(methodName), collapse=", "));
  }

  # Argument 'classNames':
  classNames <- as.character(classNames);
  if (length(classNames) == 0) {
    throw("Argument 'classNames' is empty.");
  }

  # Argument 'firstOnly':
  firstOnly <- as.logical(firstOnly);


  res <- list();
  for (kk in seq_along(classNames)) {
    className <- classNames[kk];
    fcnName <- paste(methodName, className, sep=".");
    obj <- do.call(getAnywhere, list(fcnName));
    if (length(obj$objs) == 0) {
      # No matching objects
      next;
    }

    # WORKAROUND: In R (< 3.1.?) there is a bug in getAnywhere()
    # causing it to return garbage in parts of the 'objs' list.
    hasBug <- (length(obj$objs) > length(obj$where))
    if (hasBug) {
      ## Rebuild 'objs' manually
      n <- length(obj$where)
      obj$objs <- vector("list", length=n)
      for (ii in seq_len(n)) {
        where <- obj$where[[ii]]
        tryCatch({
          if (grepl("^namespace:", where)) {
            env <- asNamespace(gsub("^namespace:", "", where))
          } else {
            env <- as.environment(where)
          }
          if (exists(fcnName, envir=env)) {
            obj$objs[[ii]] <- get(fcnName, envir=env)
          }
        }, error = function(ex) {})
      } # for (ii ...)
    }

    # Keep only functions
    keep <- which(sapply(obj$objs, FUN=is.function));
    if (length(keep) == 0) {
      # No functions
      next;
    }

    # Keep the first function
    first <- keep[1];
    fcn <- obj$objs[[first]];
    where <- obj$where[first];

    resKK <- list();
    resKK$class <- className;
    resKK$name <- methodName;
    resKK$fullname <- fcnName;
    resKK$fcn <- fcn;
    resKK$where <- obj$where;

    res[[className]] <- resKK;

    # Return only the first match?
    if (firstOnly) {
      break;
    }
  } # for (kk ...)

  res;
}, private=TRUE) # findDispatchMethodsS3()


############################################################################
# HISTORY:
# 2015-02-02
# o WORKAROUND: In R (< 3.1.?) there is a bug in getAnywhere() causing it
#   to return garbage in parts of the 'objs' list.  This bug has been
#   there all the time, but was only detected now when a package test
#   for findDispatchMethodsS3() was added.
# 2010-12-02
# o Added Rdoc comments.
# o Made findDispatchMethodsS3() a default method.
# 2009-11-20
# o Added findDispatchMethodsS3().
############################################################################
