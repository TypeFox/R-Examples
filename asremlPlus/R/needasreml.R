.onAttach <- function(libname, pkgname)
{ if(!("asreml" %in% loadedNamespaces()))
  { asreml.loaded <- requireNamespace("asreml", quietly=TRUE)
    if (!asreml.loaded)  
    { packageStartupMessage("ASReml-R needs to be loaded if the mixed-model functions are to be used.

ASReml-R is available from VSNi. Please visit http://www.vsni.co.uk/ for more information.\n")}
  }
}

