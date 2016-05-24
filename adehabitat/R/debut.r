##### Chargement de base
.onAttach <- function(libname, pkgname)
{
    msg <- paste("\n************************************************\n",
                 "************************************************\n",
                 "THE PACKAGE adehabitat IS NOW DEPRECATED!!!!!!!\n It is dangerous to use it, as bugs will no longer be corrected.\n",
                 "It is now recommended to use the packages adehabitatMA, adehabitatLT, adehabitatHR, and adehabitatHS.\n",
                 "These 4 packages are the future of adehabitat.\n They have a vignette explaining in detail how they can be used.\nThey implement more methods than adehabitat\nThey are based on the more common and more clever spatial classes implemented in sp.\nBugs are corrected frequently.\nReally, avoid to use the classical adehabitat, unless you have a very good reason for it.\n\n",
                 "*****THIS IS THE VERY LAST WARNING*****\n This is the last version of adehabitat submitted to CRAN (at the time of writing: 2015-03-27).\n THE NEXT VERSION OF adehabitat WILL JUST BE A VIRTUAL PACKAGE LOADING ALL THE PACKAGES DESCRIBED ABOVE.\n",
                 sep="")
    packageStartupMessage(msg)
}
