ifile <- function(measure, entry = 1:length(measure), ...) {
  return(as.ifile(data.frame(entry, measure, ...)))
}
