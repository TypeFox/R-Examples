.onAttach <- function(libname, pkgname){
    pkg.version <- packageDescription("bmmix", fields = "Version")
    startup.txt <- paste("\n   === bmmix", pkg.version, "is loaded ===\n")
    startup.txt <- paste(startup.txt, "Information & Documentation -> check the R-epi project: \nhttp://sites.google.com/site/therepiproject\n", sep="\n")
    startup.txt <- paste(startup.txt, "Questions -> ask the R-epi forum: \nr-epi@googlegroups.com", sep="\n")
    packageStartupMessage(startup.txt)
}
