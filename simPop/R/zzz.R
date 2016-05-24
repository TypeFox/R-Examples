.onAttach <- function(lib, pkg) {
  packageStartupMessage(paste("Package simPop",utils::packageVersion("simPop"),"has been loaded!\n"))
}