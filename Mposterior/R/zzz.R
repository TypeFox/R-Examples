.onAttach <- function (...) {
  mylib     <- dirname(system.file(package = "Mposterior"))
  ver       <- utils::packageVersion("Mposterior")
  builddate <- utils::packageDescription("Mposterior", lib.loc = mylib)$Date
  curryear  <- format(Sys.time(), "%Y")
  mess <- c("## ", "## Mposterior: R package for Robust and Scalable Bayes via a Median of Subset Posterior Measures",
            paste("## (Version ", ver,", built: ", builddate,")", sep=""),
            paste("## Copyright (C) 2013-",curryear, " Sanvesh Srivastava",sep=""))
  mess <- paste(mess, collapse = "\n")
  packageStartupMessage(mess)
}
