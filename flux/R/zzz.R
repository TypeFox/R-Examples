.onAttach <- function(lib, pkg)  {
    packageStartupMessage("This is flux ",
                          utils::packageDescription("flux", field="Version"),
                          appendLF = TRUE)
}
