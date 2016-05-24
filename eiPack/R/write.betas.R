write.betas <- function(beta, names) {
  cmd <- paste("echo \"", beta, "\" | gzip >> ", names, ".txt.gz", sep = "")
  system(cmd)
  return(invisible(0))
}
