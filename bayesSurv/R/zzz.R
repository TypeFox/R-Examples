###############################################
#### AUTHOR:    Arnost Komarek             ####
####            (2004)                     ####
####                                       ####
#### FILE:      zzz.R                      ####
####                                       ####
#### FUNCTIONS: .First.lib                 ####
###############################################

### =============================================
### .First.lib
### =============================================
.onAttach <- function(libname, pkgname)
#.First.lib <- function(lib, pkg)
{
   ###library.dynam("bayesSurv", pkgname, libname)    ## no more needed, load is provided by useDynLib in NAMESPACE

   packageStartupMessage(paste(
       "\n",
       "### Bayesian Survival Regression with Flexible Error and Random Effects Distributions \n",
       "### Arnost Komarek\n\n",
       "### See citation(\"bayesSurv\") or toBibtex(citation(\"bayesSurv\")) for the best way to cite\n",
       "### the package if you find it useful.\n\n", sep=""))   

   invisible()
}

