.onAttach <- function(lib, pkg)  {
    packageStartupMessage("This is simba ",
                          utils::packageDescription("simba", field="Version"),
                          appendLF = TRUE)
}
