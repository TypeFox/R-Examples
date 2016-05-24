.onAttach <- function(lib, pkg)  {
    packageStartupMessage("This is eHOF ",
    utils::packageDescription("eHOF", field="Version"),
    appendLF = TRUE)
    options(eHOF.bootselectmessage = TRUE)
    options(repos = c(CRAN="http://cran.r-project.org"))
}
