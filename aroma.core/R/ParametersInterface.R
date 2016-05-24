###########################################################################/**
# @RdocClass ParametersInterface
#
# @title "The ParametersInterface class interface"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/########################################################################### 
setConstructorS3("ParametersInterface", function(...) {
  extend(Interface(), "ParametersInterface");
})


###########################################################################/**
# @RdocMethod getParameters
#
# @title "Gets a list of parameters"
#
# \description{
#  @get "title" associated with the object.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a named @list.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getParameters", "ParametersInterface", function(this, ...) {
  list();
})

setMethodS3("getParameterSets", "ParametersInterface", function(this, ...) {
  paramsList <- getParameters(this, ...);

  hasSets <- attr(paramsList, "hasSets");

  # If the 'paramsList' is flat, turn it into an unnamed set of parameters
  if (is.null(hasSets)) {
    paramsList <- list(paramsList);
    hasSets <- TRUE;
    attr(paramsList, "hasSets") <- hasSets;
  }

  paramsList;  
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getParametersAsString
#
# @title "Gets the parameters as character"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getParameters".}
#  \item{collapse}{(optional) A @character string used to collapse the
#    individual parameter strings.}
#  \item{drop}{If @TRUE and there is only one set of parameters, 
#    then a single @character @vector is returned, otherwise a @list.}
# }
#
# \value{
#  Returns a @list of @character @vectors, or a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getParametersAsString", "ParametersInterface", function(this, ..., collapse=c(", ", "; "), drop=TRUE) {
  paramsList <- getParameterSets(this, ..., drop=FALSE);
  nbrOfSets <- length(paramsList);

  # If more than one set and the first one is not named and is empty,
  # then drop it.  NB: This allows subclasses to override getParameterSets()
  # instead of getParameters().
  if (nbrOfSets > 1L) {
    keys <- names(paramsList);
    # Sanity check
    if (is.null(keys)) {
      throw("INTERNAL ERROR: Detected non-named parameter sets.");
    }
    # Drop first?
    if (nchar(keys[1L]) == 0L) {
      paramsList <- paramsList[-1L];
      nbrOfSets <- length(paramsList);
    }
  }

  # Coerce each set to character strings
  res <- list();
  for (kk in seq_along(paramsList)) {
    key <- names(paramsList)[kk];
    params <- paramsList[[kk]];

    s <- trim(capture.output(str(params)))[-1L];
    s <- gsub("^[$][ ]*", "", s);
    s <- gsub(" [ ]*", " ", s);
    s <- gsub("[ ]*:", ":", s);

    if (length(collapse) > 0L) {
      s <- paste(s, collapse=collapse[1L]);
    }

    s <- sprintf("{%s}", s);

    if (!is.null(key) && nchar(key) > 0L) {
      s <- sprintf("%s=%s", key, s);
    }

    res[[kk]] <- s;
  } # for (kk ...)

  if (drop && nbrOfSets == 1L) {
    res <- res[[1L]];
  } else if (nbrOfSets > 1L && length(collapse) > 1L) {
    res <- paste(res, collapse=collapse[2L]);
  }

  res;
}) # getParametersAsString()


############################################################################
# HISTORY:
# 2012-11-21
# o Now getParametersAsString() handles sets of parameters as well.
# o Added getParameterSets() to ParametersInterface.
# 2012-11-20
# o Added Rdoc comments.
# o Created.
############################################################################
