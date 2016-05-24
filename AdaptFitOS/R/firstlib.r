.onAttach <- function (lib, pkg) {
  ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"),"Version")
  ver <- as.character(ver)
  packageStartupMessage("\nAdaptFitOS ",ver," loaded. Type 'help(\"AdaptFitOS-package\")' for an overview.\n",
                        "\nPlease cite as:\n   Wiesenfarth, M., Krivobokova, T., Klasen, S., & Sperlich, S. (2012).\n   Direct simultaneous inference in additive models and its application to model undernutrition.\n   Journal of the American Statistical Association, 107(500), 1286-1296."                     
                         , domain = NULL,  appendLF = TRUE)
}

.onLoad <- function(...) {
}
