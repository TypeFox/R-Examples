###########################################################################/**
# @RdocClass ColumnNamesInterface
#
# @title "The ColumnNamesInterface class interface"
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
setConstructorS3("ColumnNamesInterface", function(...) {
  extend(Interface(), "ColumnNamesInterface");
})



###########################################################################/**
# @RdocMethod nbrOfColumns
#
# @title "Gets the number of columns"
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
#   Returns an @integer.
#   If the number of columns cannot be inferred, @see NA is returned.
# }
# @author
#
# \seealso{
#   @seemethod "getColumnNames".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("nbrOfColumns", "ColumnNamesInterface", function(this, ...) {
  columns <- getColumnNames(this);
  if (is.null(columns)) return(NA_integer_);
  length(columns);
})



###########################################################################/**
# @RdocMethod getDefaultColumnNames
#
# @title "Gets the default column names"
#
# \description{
#   @get "title", that is, the column names without translations.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "getColumnNames".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getDefaultColumnNames", "ColumnNamesInterface", abstract=TRUE, protected=TRUE);



###########################################################################/**
# @RdocMethod getColumnNames
#
# @title "Gets the column names"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{translate}{If @TRUE and a names translator is set, the
#     column names are translated before returned.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "getDefaultColumnNames".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getColumnNames", "ColumnNamesInterface", function(this, ..., translate=TRUE) {
  names <- getDefaultColumnNames(this, ...);

  # Translate?
  if (translate) {
    names <- translateColumnNames(this, names);
  }

  names;
})



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TRANSLATOR FUNCTIONS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("clearListOfColumnNamesTranslators", "ColumnNamesInterface", function(this, ...) {
  setListOfColumnNamesTranslators(this, list());
}, protected=TRUE)

setMethodS3("clearColumnNamesTranslator", "ColumnNamesInterface", function(this, ...) {
  clearListOfColumnNamesTranslators(this);
})


setMethodS3("getListOfColumnNamesTranslators", "ColumnNamesInterface", function(this, ...) {
  res <- this$.listOfColumnNamesTranslators;
  if (is.null(res)) {
    res <- list();
  }
  res;
}, protected=TRUE)

setMethodS3("setListOfColumnNamesTranslators", "ColumnNamesInterface", function(this, fnList, ...) {
  # Argument 'fnList':
  for (kk in seq_along(fnList)) {
    fcn <- fnList[[kk]];
    if (!is.function(fcn)) {
      throw("Element #", kk, " of argument 'fnList' is not a function: ",
                                                           class(fcn)[1]);
    }
  }

  this$.listOfColumnNamesTranslators <- fnList;

  invisible(this);
}, protected=TRUE)


setMethodS3("getColumnNamesTranslator", "ColumnNamesInterface", function(this, ...) {
  fnList <- getListOfColumnNamesTranslators(this, ...);

  # No names translator?
  if (length(fnList) == 0) {
    return(NULL);
  }

  # Create names translator function
  res <- function(names, ...) {
    for (kk in seq_along(fnList)) {
      fcn <- fnList[[kk]];
      names <- fcn(names, ...);
    }
    names;
  }
  res;
}, protected=TRUE)



setMethodS3("translateColumnNames", "ColumnNamesInterface", function(this, names, ...) {
  nameTranslator <- getColumnNamesTranslator(this);

  if (!is.null(nameTranslator)) {
    names2 <- nameTranslator(names, file=this);

    # Sanity check
    if (any(is.na(names2))) {
      throw("Failed to translate names. Some names were translated to NA:s ",
            paste(head(names[is.na(names2)]), collapse=", "));
    }
    if (length(names2) != length(names)) {
      throw(sprintf("Failed to translate column names. The translator is erroneous, because it drops/adds some names (passed %d names but got %d names).", length(names), length(names2)));
    }
    names <- names2;

    if (identical(attr(names, "isFinal"), TRUE))
      return(names);
  }

  # Do nothing
  names;
}, private=TRUE)


setMethodS3("appendColumnNamesTranslatorByNULL", "ColumnNamesInterface", function(this, ...) {
  # Nothing to append
  invisible(this);
}, protected=TRUE)


setMethodS3("appendColumnNamesTranslatorBylist", "ColumnNamesInterface", function(this, list, ...) {
  # Arguments 'list':
  if (!inherits(list, "list")) {
    throw("Argument 'list' is not a list: ", class(list)[1]);
  }

  for (kk in seq_along(list)) {
    by <- list[[kk]];
    appendColumnNamesTranslator(this, by, ...);
  }
}, protected=TRUE)


setMethodS3("appendColumnNamesTranslatorBycharacter", "ColumnNamesInterface", function(this, names, ...) {
  # Validate argument 'names'
  names <- Arguments$getCharacters(names);

  # Append a translator function that always returns a constant string
  appendColumnNamesTranslator(this, function(...) { names });
}, protected=TRUE)


setMethodS3("appendColumnNamesTranslatorByfunction", "ColumnNamesInterface", function(this, fcn, ..., validate=TRUE) {
  # Arguments 'fcn':
  if (!is.function(fcn)) {
    throw("Argument 'fcn' is not a function: ", class(fcn)[1]);
  }

  # Sanity check
  if (validate) {
    names <- getDefaultColumnNames(this);
    namesT <- fcn(names, file=this);

    # More sanity checks
    if (length(namesT) != length(names)) {
      throw(sprintf("Argument 'fcn' specifies a translator function that return %d string(s) when give %d: %s", length(names), length(namesT), hpaste("'", namesT, "'")));
    }
  }

  fnList <- getListOfColumnNamesTranslators(this);
  fnList <- c(fnList, fcn);
  setListOfColumnNamesTranslators(this, fnList);
}, protected=TRUE)


setMethodS3("appendColumnNamesTranslator", "ColumnNamesInterface", function(this, by, ...) {
  # Arguments 'by':
  classNames <- class(by);
  methodNames <- sprintf("appendColumnNamesTranslatorBy%s", classNames);

  keep <- sapply(methodNames, FUN=exists, mode="function");
  methodNames <- methodNames[keep];

  if (length(methodNames) == 0) {
    throw("Failed to set the names translator. Could not find an appendColumnNamesTranslatorBy<className>() function for this object: ", paste(classNames, collapse=", "));
  }

  methodName <- methodNames[1];
  fcn <- get(methodName, mode="function");
  res <- fcn(this, by, ...);

  # Allow the object to update itself according to these new rules.
  updateColumnNames(this);

  invisible(res);
}, protected=TRUE)


setMethodS3("setColumnNamesTranslator", "ColumnNamesInterface", function(this, ...) {
  clearListOfColumnNamesTranslators(this);
  appendColumnNamesTranslator(this, ...);
})



###########################################################################/**
# @RdocMethod setColumnNames
#
# @title "Sets the column names"
#
# \description{
#   @get "title".
#   This is done using a names translator function that returns the
#   specified names.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments, typically a @character string, which are
#     passed to the names translator generator.
#  }
# }
#
# \value{
#   Returns (invisibly) itself.
# }
#
# @author
#
# \seealso{
#   @seemethod "getColumnNames".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setColumnNames", "ColumnNamesInterface", function(this, ...) {
  # Set a translator function that always returns a constant
  setColumnNamesTranslator(this, ...);
})



setMethodS3("updateColumnNames", "ColumnNamesInterface", function(this, ...) {
}, protected=TRUE)




############################################################################
# HISTORY:
# 2013-12-18
# o Now nbrOfColumns() for ColumnNamesInterface returns NA if column names
#   cannot be inferred and hence not be counted.
# 2012-11-07
# o Now clearListOfColumnNamesTranslators() utilizes
#   setListOfColumnNamesTranslators().
# 2012-11-02
# o Created from FullNameInterface.R
############################################################################
