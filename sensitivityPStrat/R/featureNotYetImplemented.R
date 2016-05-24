.FeatureNotYetImplemented <- function(feature)
  stop(gettextf("'%s' is not implemented yet", paste(as.character(sys.call(sys.parent())[[1L]]))), call. = FALSE)
