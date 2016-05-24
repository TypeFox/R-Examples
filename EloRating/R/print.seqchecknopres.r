print.seqchecknopres <- function(x, ...) {
  if(x$checksum["presence"]==1) cat("no presence data supplied\n")
  
  if(x$checksum["IDcheck"]==1) cat("IDs occur in the data with inconsistent capitalization (ignore if on purpose)...WARNING\n")
  if(x$checksum["selfinteractions"]==1) cat(x$selfinteractions, "\n", sep="")
  if(x$checksum["length"]==1) cat("your data vectors do not match in length...ERROR\n")
  if(x$checksum["singledayobs"]==1) cat("the following individuals were observed only on one day: ", paste(x$singledaycases, collapse=", "), " ...WARNING\n", sep="" )
  
  if(x$checksum["selfinteractions"]+x$checksum["IDcheck"]+x$checksum["length"]+x$checksum["singledayobs"]==0) cat("Everything seems to be fine with the interaction sequence...OK\n")
}