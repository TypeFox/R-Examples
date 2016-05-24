touch.betas <- function(names) {
  check <- system("ls", intern = TRUE)
  if (any(check %in% names)) {
    exist <- na.omit(check[match(names, check)])
    stop(paste( length(exist), " files already exist in ", getwd(), ".  Either select a new working directory or remove these files and try again.", sep = ""))
  }
  for (ii in names) {
    cmd <- paste("echo \"", ii, "\" | gzip >> ", ii, sep = "")
    system(cmd)
  }
  return(invisible(0))
}
