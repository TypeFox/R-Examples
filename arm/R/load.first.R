.onAttach <- function(...) {
  mylib <- dirname(system.file(package = "arm"))
  ver <- packageDescription("arm", lib.loc = mylib)$Version
  builddate <- packageDescription("arm", lib.loc = mylib)$Date
  packageStartupMessage(paste("\narm (Version ", ver, ", built: ", builddate, ")\n", sep = ""))
  packageStartupMessage("Working directory is ", getwd(), "\n")
}
