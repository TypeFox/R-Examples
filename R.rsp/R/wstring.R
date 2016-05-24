setMethodS3("wstring", "default", function(..., sep="", envir=parent.frame()) {
  s <- paste(..., sep=sep);

  # Nothing to do?
  if (length(s) == 0L) return(s);
  if (length(s) > 1L) {
    Recall <- sys.function();  # base::Recall() does not work with *apply()
    return(sapply(s, FUN=Recall, envir=envir));
  }

  # Nothing to do?
  if (regexpr("{{", s, fixed=TRUE) == -1L) {
    return(s);
  }

  bfr <- NULL;
  pattern <- "{{(.*?)}}";
  while ((pos <- regexpr(pattern, s, perl=TRUE)) != -1L) {
    # Parse
    len <- attr(pos, "match.length");
    head <- substring(s, first=1L, last=pos-1L);
    code <- substring(s, first=pos+2L, last=pos+len-3L);
    s <- substring(s, first=pos+len);

    # Evaluate
    expr <- parse(text=code);
    value <- eval(expr, envir=envir);
    value <- as.character(value);

    bfr <- c(bfr, head, value);
  }
  bfr <- c(bfr, s);
  bfr <- paste(bfr, collapse="");

  bfr;
}, protected=TRUE) # wstring()


############################################################################
# HISTORY:
# 2013-12-21
# o Added wstring().  Should probably endup in R.utils one day.
# o Created.
############################################################################
