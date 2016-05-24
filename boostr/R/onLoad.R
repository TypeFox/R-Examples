# .First <- function() {
#   options(
#     prompt="> ",
#     repos = c(CRAN = "http://cran.rstudio.com/"),
#     browserNLdisabled = TRUE,
#     deparse.max.lines = 2)
#   
#   if (interactive()) {
#     # suppressMessages(require(devtools))
#   }
#   
#   # Sys.setenv(PATH=paste0(Sys.getenv("PATH"), ":/home/steven/.cabal/bin"))
# }
# 
# 
# .onLoad <- function(libname, pkgname) {
# #   .First()
#   rmdWeave <- function(file, ...) {
#     rmarkdown::render(file, "all")
#   }
#   tools::vignetteEngine("rmarkdown", weave = rmdWeave,
#                         tangle = function(file, ...) {file},
#                         pattern = "[.]Rmd$", package = "rmarkdown")
# }