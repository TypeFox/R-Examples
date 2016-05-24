.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste("synchrony", 
          utils::packageDescription("synchrony",
                                    field="Version"),
          "loaded."),
    appendLF = TRUE)
}