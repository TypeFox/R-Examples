.onAttach <- function (lib, pkg){
  packageStartupMessage("This is biotic, version ",
                        utils::packageDescription("biotic", fields="Version"),
                        appendLF = TRUE)
}
