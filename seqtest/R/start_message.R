#' @importFrom utils packageDescription
#'
.onAttach <- function(libname, pkgname){


  desc <- packageDescription("seqtest")
  d1 <- desc$Version
  nk <- paste0(rep(" ", 6 - nchar(d1)))

  packageStartupMessage("|----------------------------|\n",
                          paste0("| ", desc$Package, " ", d1," (",desc$Date,")"), nk, "|\n" ,
                        "| Sequential Triangular Test |\n" ,
                        "|----------------------------|" )

}

version <- function(pkg = "seqtest") {

  lib <- dirname(system.file(package = pkg))
  desc <- packageDescription(pkg)

  return(paste(desc$Package, desc$Version, desc$Date,lib))

}
