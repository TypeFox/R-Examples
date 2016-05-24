pfile <- function(measure, entry = 1:length(measure), ...) {
  return(as.pfile(data.frame(entry, measure, ...)))
}
