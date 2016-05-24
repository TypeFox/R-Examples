.onLoad <- function(lib, pkg){
#    require("methods", character = TRUE, quietly = TRUE)
}


.onAttach <- function(library, pkg)
{
buildStartupMessage(pkg = "distrTeach", library = library, packageHelp = TRUE, 
#                    MANUAL="http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf",
VIGNETTE = gettext("Package \"distrDoc\" provides a vignette to this package as well as to several related packages; try vignette(\"distr\")."))
  invisible()
}



