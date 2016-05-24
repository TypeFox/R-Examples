library("R.methodsS3")
library("R.oo")
library("R.utils")

setConstructorS3("ROptions", function(...) {
  extend(Options(), "ROptions",
    .helpText = NULL,
    .descriptions = list()
  );
})

setMethodS3("as.list", "ROptions", function(this, ...) {
  # Import R options each time.
  this$.options <- options();
  tryCatch({
    NextMethod("as.list");
  }, error = function(ex) {
    # In case this method was called explicitly, e.g. ROptions$as.list()
    as.list.Options(this);
  })
})

setMethodS3("setOption", "ROptions", function(this, ...) {
  oldValue <- NextMethod("setOption");
  # Store options
  options(this$.options);
  invisible(oldValue);
})


setMethodS3("getDescription", "ROptions", function(this, option, default="", force=FALSE, ...) {
  # Is option value cached?
  if (!force) {
    value <- ROptions$.descriptions[[option]];
    if (!is.null(value))
      return(value);
  }

  # Search for option value in help(options)
  helpText <- this$.helpText;
  if (force || is.null(helpText)) {
    helpText <- readRdHelp("options");
    helpText <- trim(gsub("_\b", "", helpText));
    this$.helpText <- helpText;
  } else {
    helpText <- this$.helpText;
  }

  # Pattern to search for
  match <- paste("'", option, "':", sep="");

  # Find start position
  pos <- regexpr(match, helpText);
  idx <- which(pos != -1);
  if (length(idx) == 0)
    return(default);
  idx <- idx[1];
  value <- helpText[idx:length(helpText)];

  # Find stop position
  idx <- which(nchar(value) == 0)[1];
  value <- value[1:(idx-1)];
  value <- paste(value, collapse=" ");
  value <- sub(match, "", value);
  value <- trim(value);

  # Store value
  ROptions$.descriptions[[option]] <- value;

  value;
})



setMethodS3("setDescription", "ROptions", function(this, option=NULL, value=NULL, collapse="", sep="", ...) {
  oldValue <- getDescription(this, option);

  value <- as.character(value);
  value <- paste(value, collapse=collapse, sep=sep);
  descriptions <- ROptions$.descriptions;
  ROptions$.descriptions[[option]] <- value;

  invisible(oldValue);
})



setMethodS3("getDataTypes", "ROptions", function(this, option, default="", force=FALSE, ...) {
  # Is option value cached?
  if (!force) {
    value <- ROptions$.dataTypes[[option]];
    if (!is.null(value))
      return(value);
  }

  # Check current value and return that as the accepted mode.
  optValue <- getOption(this, option);
  if (!is.null(optValue))
    return(storage.mode(optValue));

  # Finally, if nothing is found, return NULL.
  NULL;
})



rOptions <- ROptions();

###########################################################################
# HISTORY:
# 2013-09-18
# o Code no longer assumes that packages R.methodsS3, R.oo and R.utils
#   are attached.
###########################################################################
