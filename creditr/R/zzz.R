.onLoad <- function(libname, pkgname) {
  
  ## compress the vignette PDF to fix CMD Check WARNING
  ## The below works and please don't delete this!
  
  Sys.setenv("_R_BUILD_COMPACT_VIGNETTES_" = "qpdf")
}

.onAttach <- function(...) {
  packageStartupMessage(
    "\nIF YOU USE THIS CDS PACKAGE, YOUR USE WILL SIGNIFY YOUR
    UNDERSTANDING \nAND IRREVOCABLE ACCEPTANCE OF THIS LICENSE AND ITS
    TERMS, WHICH INCORPORATE \nBY REFERENCE THE INTERNATIONAL SWAPS AND
    DERIVATIVES ASSOCIATION, INC.'S CDS \nSTANDARD MODEL PUBLIC LICENSE,
    WHICH IS AVAILABLE AT
    \nhttp://www.cdsmodel.com/cdsmodel/cds-disclaimer.page before using
    the package. \n\nNOTHING IN THIS LICENSE RESTRICTS YOUR ABILITY TO USE
    THE ISDA(R) CDS \nSTANDARD MODEL.\n \nDISCLAIMER: ISDA HAS NEITHER
    REVIEWED, APPROVED NOR ENDORSED THE USE OF \nTHE CDS PACKAGE. THOSE
    PERSONS USING THIS CDS PACKAGE ARE ENCOURAGED TO SEEK \nTHE ADVICE OF
    A LEGAL PROFESSIONAL TO EVALUATE ITS SUITABILITY FOR THEIR USE.\n
    \nISDA(R) is a registered mark of the International Swaps and
    Derivatives \nAssociation, Inc.\n \nPlease type **yes** to assent to
    the aforementioned terms.\n")
  
  if(interactive()){
    while (readLines(n=1)!="yes")
      packageStartupMessage("Please type **yes** to assent to the
                             aforementioned terms.\n")
  }
}
