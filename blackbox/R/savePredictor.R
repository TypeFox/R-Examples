# probably completely obsolete since predict.OKrig -> Cpredict -> objects in memory of dll
savePredictor <- function(file) {
  if (missing(file)) {
    stop.redef("'file' argument missing, with no default")
  }
  unlink(file) ## this *deletes* the file if it exists
  save(  predictor = blackbox.getOption("fitobject") , file=file)
  msgg <- paste("\"", file, "\"", sep="")
  msg <- paste("Use load(", msgg, ") to load 'predictor' in a future R session", sep="")
  cat(msg) ## with print() all quotation mark characters are expressed as "
  invisible(NULL)
}
