#*** zzz.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
#* ********************************************************************************* */

.onAttach <- function(libname, pkgname)
#.First.lib <- function(libname, pkgname)
{
   ##library.dynam("mixAK", pkgname, libname)   ## no more needed, load is provided by useDynLib in the NAMESPACE

   packageStartupMessage(paste(
       "\n",
       "### Mixture of methods including mixtures\n",
       "### Arnost Komarek\n\n",
       "### See citation(\"mixAK\") or toBibtex(citation(\"mixAK\")) for the best way to cite\n",
       "### the package if you find it useful.\n\n", sep=""))
   #cat("\n")
   #cat("### Mixture of methods including mixtures\n")
   #cat("### Arnost Komarek\n\n")
   #cat("### See citation(\"mixAK\") or toBibtex(citation(\"mixAK\")) for the best way to cite\n")
   #cat("### the package if you find it useful.\n\n")
   
   invisible()
}

