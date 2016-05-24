ConfigShowDefinitions <- function(configuration) {
  defaults <- grep("^DEFAULT_", names(configuration$definitions))
  invisible(lapply(as.vector(paste(names(configuration$definitions)[defaults], paste(unlist(configuration$definitions)[defaults], "\n", sep = ''), sep = " = ")), cat))
  invisible(lapply(as.vector(paste(names(configuration$definitions)[-defaults], paste(unlist(configuration$definitions)[-defaults], "\n", sep = ''), sep = " = ")), cat))
}
