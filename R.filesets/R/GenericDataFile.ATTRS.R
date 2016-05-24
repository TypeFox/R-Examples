setMethodS3("getAttributes", "GenericDataFile", function(this, ...) {
  attrs <- this$.attributes;
  if (length(attrs) == 0) {
    attrs <- list();
  } else {
    # Always return attributes in lexicographic order by names
    names <- names(attrs);
    if (length(names) > 0) {
      o <- order(names);
      attrs <- attrs[o];
    }
  }
  attrs;
}, protected=TRUE)



setMethodS3("setAttributes", "GenericDataFile", function(this, ...) {
  # Argument '...':
  args <- list(...);
  names <- names(args);
  if (is.null(names)) {
    throw("No named arguments specified.");
  }
  
  # Update the attributes.
  attrs <- this$.attributes;
  attrs[names] <- args;
  this$.attributes <- attrs;

  invisible(args);
}, protected=TRUE)



setMethodS3("getAttribute", "GenericDataFile", function(this, name, defaultValue=NULL, ...) {
  attrs <- this$.attributes;
  if (name %in% names(attrs)) {
    value <- attrs[[name]];
  } else {
    value <- defaultValue;
  }
  value;
}, protected=TRUE)



setMethodS3("setAttribute", "GenericDataFile", function(this, name, value, ...) {
  attrs <- this$.attributes;
  attrs[[name]] <- value;
  this$.attributes <- attrs;

  invisible(attrs[name]);
}, protected=TRUE)



setMethodS3("testAttributes", "GenericDataFile", function(this, select, ...) {
  # Get the attributes to be tested
  attrs <- getAttributes(this);
  expr <- substitute(select);
  res <- eval(expr, envir=attrs, enclos=parent.frame());
  res;
}, protected=TRUE)



setMethodS3("setAttributesBy", "GenericDataFile", function(this, object, ...) {
  if (inherits(object, "character")) {
    setAttributesByTags(this, object, ...);
  } else {
    throw("Unknown type on argument 'object': ", class(object)[1]);
  }
}, protected=TRUE)



setMethodS3("setAttributesByTags", "GenericDataFile", function(this, tags=getTags(this), ...) {
  # Split tags
  if (length(tags) > 0) {
    tags <- unlist(strsplit(tags, split=","), use.names=FALSE);
    tags <- trim(tags);
  }

  newAttrs <- list();

  # Get all <name>=<value> tags
  pattern <- "^([abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ]+)=(.*)$";
  values <- grep(pattern, tags, value=TRUE);
  for (kk in seq_along(values)) {
    tag <- values[[kk]];
    key <- gsub(pattern, "\\1", tag);
    value <- gsub(pattern, "\\2", tag);

    # Try to coerce:
    suppressWarnings({
      value2 <- as.integer(value);
      if (!identical(value2 == value, TRUE)) {
        value2 <- as.double(value);
        if (!identical(value2 == value, TRUE)) {
          value2 <- as.character(value);
        }
      }
      value <- value2;
    })

    newAttrs <- c(newAttrs, setAttribute(this, key, value));
  }

  # Return updated attributes
  invisible(newAttrs);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2008-12-10
# o BUG FIX: getAttributes() for GenericDataFile:s would give an error
#   if there were no attributes.
# 2007-05-09
# o Moved all attribute features from AffymetrixFile to GenericDataFile.
# 2007-03-06
# o Added setAttributesBy().
# 2007-03-05
# o Added setAttributesByTags(), which now also tries to coerce values.
# o Added support for (in-memory) attributes.
# 2007-02-07
# o Added getChecksum(), writeChecksum(), readChecksum(), and 
#   compareChecksum() and validateChecksum(). I did this because I noticed 
#   by chance that some of my CEL files transferred via an external HDD got
#   corrupt probe signals.
# 2007-01-14
# o Added a test for "unknown" (=unused) arguments to constructor.
# 2007-01-07
# o Added hasTags() and hasTag().
# 2006-11-02
# o Added getFullName(), getTags() and redefined getName().
# 2006-09-15
# o Added stextSize().
# 2006-08-27
# o Added stextLabel() and stextLabels(). stext is for "side text", cf. 
#   mtext for "margin text". stext() is slightly more convenient than mtext
#   when it comes to different font sizes.
# o Added copyTo().
# 2006-08-14
# o Added abstract fromFile().
# 2006-08-11
# o Created from AffymetrixDataFile in order to represent CDF files too.
############################################################################
