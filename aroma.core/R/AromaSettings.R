setConstructorS3("AromaSettings", function(basename=NULL, ...) {
  if (inherits(basename, "Settings")) {
    this <- basename;
    class(this) <- c("AromaSettings", class(this));
  } else {
    this <- extend(Settings(basename=basename, ...), "AromaSettings");
  }
  this;
})


setMethodS3("getVerbose", "AromaSettings", function(this, default=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Transition rule 'argVerbose':
  # Methods have argument verbose=FALSE by default.  If not overridden,
  # then argument 'default' of this method will be identical to ditto.
  # When this happens, we can choose to either:
  #  (i) keep interpreting this as is (old style), or
  # (ii) have the value to default to the settings.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  transRule <- getOption(this, "transitionRules/useVerboseOption", FALSE);
  if (transRule) {
    callDepth <- sys.nframe()-1L;
    maxCallDepth <- getOption(this, "verbose/maxCallDepth", 15);
    if (callDepth <= maxCallDepth) {
      verbose <- getOption(this, "verbose", default);
    } else {
      verbose <- default;
    }
  } else {
    verbose <- default;
  }

  verbose <- Arguments$getVerbose(verbose);
  verbose;
})


setMethodS3("setVerbose", "AromaSettings", function(this, ..., timestamp=TRUE) {
  verbose <- Arguments$getVerbose(..., timestamp=timestamp);
  setOption(this, "verbose", verbose);
})


setMethodS3("getRam", "AromaSettings", function(this, default=1, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Transition rule 'argRam':
  # Methods have argument ram=NULL (or ram=1) by default.  If not overridden,
  # then argument 'default' of this method will be identical to ditto.
  # When this happens, we can choose to either:
  #  (i) keep interpreting this as is (old style), or
  # (ii) have the value to default to the settings.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(default)) {
    default <- 1;
  }
  transRule <- getOption(this, "transitionRules/useRamOption", TRUE);
  if (transRule) {
    ram <- getOption(this, "memory/ram", default);
  } else {
    ram <- default;
  }

  ram <- Arguments$getDouble(ram, range=c(0.001, Inf));

  ram;
})


setMethodS3("setRam", "AromaSettings", function(this, value=1, ...) {
  # Argument 'value':
  value <- Arguments$getDouble(value, range=c(0.001, Inf));

  setOption(this, "memory/ram", value);
})



############################################################################
# HISTORY:
# 2012-11-26
# o BUG FIX: getRam() and setRam() for AromaSettings did not use
#   'memory/ram'.
# 2009-02-22
# o Added argument {get|set}Ram().
# o Added argument {get|set}Verbose().
# o Created.
############################################################################
