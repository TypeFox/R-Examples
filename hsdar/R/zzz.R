.onAttach <- function(lib, pkg)  {
    packageStartupMessage("\n###################################",
                          "\n  This is hsdar ",
                          utils::packageDescription("hsdar",
                                                    field="Version"),
                          "\n  To get citation entry type",
                          "\n      'citation(\"hsdar\")'",
                          "\n###################################",
                          appendLF = TRUE)
}
