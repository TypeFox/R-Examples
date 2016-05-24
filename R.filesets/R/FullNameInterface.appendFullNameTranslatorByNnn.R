setMethodS3("appendFullNameTranslatorBycharacter", "FullNameInterface", function(this, fullname, ...) {
  # Validate argument 'fullname'
  fullname <- Arguments$getCharacter(fullname, length=c(1,1));

  # Append a translator function that always returns a constant string
  appendFullNameTranslator(this, function(...) { fullname });
}, protected=TRUE)


setMethodS3("appendFullNameTranslatorByfunction", "FullNameInterface", function(this, fcn, ...) {
  # Arguments 'fcn':
  if (!is.function(fcn)) {
    throw("Argument 'fcn' is not a function: ", class(fcn)[1]);
  }

  # Sanity check
  name <- c("foo bar");
  nameT <- fcn(name, file=this);

  # More sanity checks
  if (length(nameT) != 1L) {
    throw("Argument 'fcn' specifies a translator function that does not return exactly one string if given one string: ", length(nameT));
  }

  fnList <- getListOfFullNameTranslators(this);
  fnList <- c(fnList, fcn);
  setListOfFullNameTranslators(this, fnList);
}, protected=TRUE)


setMethodS3("appendFullNameTranslatorBydata.frame", "FullNameInterface", function(this, df, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'df':
  if (!is.data.frame(df)) {
    throw("Argument 'df' is not a data.frame: ", class(df)[1]);
  }

  colnames <- colnames(df);

  reqColNamesList <- list(
    fixed=c("fixed", "replacement"),
    pattern=c("pattern", "replacement")
  );

  if (is.null(colnames) && ncol(df) == 2) {
    colnames <- reqColNamesList[["pattern"]];  # Assume pattern
    colnames(df) <- colnames;
  }

  keep <- sapply(reqColNamesList, FUN=function(x) {
    all(is.element(x, colnames));
  });
  keep <- which(keep);

  if (length(keep) == 0) {
    d <- sapply(reqColNamesList, FUN=function(s) {
      paste(sprintf("'%s'", s), collapse=", ");
    });
    d <- sprintf("(%s)", d);
    msg <- sprintf("The specified data frame does not have all of the required columns (%s): %s", paste(d, collapse=" OR "), paste(colnames, collapse=", "));
    throw(msg);
  }

  flavor <- names(keep);
  reqColNames <- reqColNamesList[[flavor]];

  lookup <- reqColNames[1];

  if (flavor == "fixed") {
    fixed <- TRUE;
  } else if (flavor == "pattern") {
    fixed <- FALSE;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate regular expression
  patterns <- df[,lookup];
  patterns <- as.character(patterns);
  replacements <- df[,"replacement"];
  replacements <- as.character(replacements);
  nbrOfRules <- length(patterns);

  # Generate translator function
  fcn <- function(names, ...) {
    # For each rule
    for (kk in seq_len(nbrOfRules)) {
      pattern <- patterns[kk];
      idxs <- grep(pattern, names, fixed=fixed);
      # No matches?
      if (length(idxs) == 0)
        next;

      # Translate
      replacement <- replacements[kk];
      names[idxs] <- gsub(pattern, replacement, names[idxs], fixed=fixed);
    } # for (kk ...)

    # Drop empty tags
    names <- gsub("[,]+", ",", names, fixed=FALSE);
    names <- gsub(",$", "", names, fixed=FALSE);

    names;
  } # fcn()

  appendFullNameTranslator(this, fcn);
}, protected=TRUE)


setMethodS3("appendFullNameTranslatorByTabularTextFile", "FullNameInterface", function(this, df, ...) {
  # Arguments 'df':
  if (!inherits(df, "TabularTextFile")) {
    throw("Argument 'df' is not a TabularTextFile: ", class(df)[1]);
  }

  df <- readDataFrame(df, defColClass="character");

  appendFullNameTranslator(this, df, ...);
}, protected=TRUE)


setMethodS3("appendFullNameTranslatorByTabularTextFileSet", "FullNameInterface", function(this, ds, ...) {
  # Arguments 'ds':
  if (!inherits(ds, "TabularTextFileSet")) {
    throw("Argument 'ds' is not a TabularTextFileSet: ", class(ds)[1]);
  }

  dummy <- sapply(ds, function(df) {
    appendFullNameTranslator(this, df, ...);
  });

  invisible(this);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2011-03-11
# o Now appendFullNameTranslatorBy<what>() for <character> and <function>
#   assert that the translator correctly returns exactly one string.
# 2010-10-17
# o Now appendFullNameTranslator(..., df) for FullNameInterface takes
#   either 'pattern' or 'fixed' translations in data.frame.
# 2010-05-26
# o Added appendFullNameTranslatorBy...() method for TabularTextFileSet:s.
# 2010-05-25
# o Added appendFullNameTranslatorBy...() method for data frames and
#   TabularTextFile:s.
# o Moved to its own file.
############################################################################
