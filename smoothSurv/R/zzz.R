###############################################
#### AUTHOR:    Arnost Komarek             ####
####            (2003)                     ####
####                                       ####
#### FILE:      zzz.R                      ####
####                                       ####
#### FUNCTIONS: .First.lib                 ####
###############################################

### =============================================
### .First.lib
### =============================================
.onAttach <- function(libname, pkgname)
#.First.lib <- function(libname, pkgname)
{
   ###library.dynam("smoothSurv", pkgname, libname)    ### no more needed, load is provided by useDynLib in NAMESPACE

   packageStartupMessage(paste(
       "\n",
       "### Survival Regression with Smoothed Error Distribution \n",
       "### Arnost Komarek\n\n",
       "### See citation(\"smoothSurv\") or toBibtex(citation(\"smoothSurv\")) for the best way to cite\n",
       "### the package if you find it useful.\n\n", sep=""))
   
   invisible()
}
