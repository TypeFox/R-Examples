.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
      packageStartupMessage("It is highly recommended that WGCNA be installed.")
      packageStartupMessage("To do this, enter the following commands into the console:")
      packageStartupMessage("setRepositories(ind=1:2)")
      packageStartupMessage("install.packages(\"WGCNA\")")
      packageStartupMessage("source(\"http://bioconductor.org/biocLite.R\")")
      packageStartupMessage("biocLite(\"AnnotationDbi\", type=\"source\")")
      packageStartupMessage("biocLite(\"GO.db\")")
      packageStartupMessage("If this fails see http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/ \n for further information.")
    }
}


