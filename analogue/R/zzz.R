.onAttach <- function(lib, pkg)  {
    packageStartupMessage("analogue version ",
                          utils::packageDescription("analogue",
                                                    fields ="Version"),
                          appendLF = TRUE)
}
