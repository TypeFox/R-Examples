.textPrompt <- function(prompt, options=c("Y"="yes", "n"="no"), caseSensitive=FALSE, maxTries=Inf, type=c("message", "output"), onSink=c("error", "ignore"), ...) {
  # Argument 'maxTries':
  maxTries <- Arguments$getNumeric(maxTries, range=c(1,Inf));

  # Argument 'onSink':
  onSink <- match.arg(onSink);

  # Argument 'type':
  type <- match.arg(type);

  # Check if standard output is redirected.
  # NOTE: There always *2* sinks for type="message", cf. help("sink").
  minSink <- switch(type, message=2L, 0L);
  hasSink <- (sink.number(type=type) > minSink);
  if (hasSink) {
    if (onSink == "error") {
      throw("Cannot prompt user via the standard ", type, ", because it is currently redirected and (most likely) not visible to the user.");
    }
  }

  # How to present the options to the user
  keys <- names(options);
  if (is.null(keys)) {
    keys <- options;
  }
  promptF <- sprintf("%s [%s]: ", prompt, paste(keys, collapse="/"));
  promptR <- gettext("Unknown reply.\n");

  # Where to prompt
  if (type == "output") {
    con <- stdout();
  } else {
    con <- stderr();
  }

  count <- 0L;
  while (count < maxTries) {
    cat(file=con, promptF);
    ans <- readline();

    count <- count + 1L;

    ans <- trim(ans);
    if (ans == "") {
      idx <- 1L;
    } else {
      if (caseSensitive) {
        idx <- pmatch(ans, options);
      } else {
        idx <- pmatch(tolower(ans), tolower(options));
      }
    }

    if (is.finite(idx)) {
      ans <- options[idx];
      break;
    }

    cat(file=con, promptR);
  } # while(...)

  attributes(ans) <- NULL;

  ans;
} # .textPrompt()


############################################################################
# HISTORY:
# 2014-05-01
# o Now argument 'type' defaults to c("message", "output").
# 2011-12-30
# o Added textPrompt().
############################################################################
