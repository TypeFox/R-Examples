.onAttach <- function(lib, pkg) {
  packageStartupMessage(paste("Package sdcTable",utils::packageVersion("sdcTable"),"has been loaded!\n"))
}
