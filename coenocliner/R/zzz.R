.onAttach <- function(lib, pkg)  {
    packageStartupMessage("This is coenocliner ",
                          utils::packageDescription("coenocliner",
                                                    field = "Version"),
                          appendLF = TRUE)
}
