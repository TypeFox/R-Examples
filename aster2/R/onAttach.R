
.onAttach <- function(lib, pkg) {
     packageStartupMessage(
         "This is beta software.\n",
         "Unless you need to do aster models with dependence groups,\n",
         "    use package \"aster\" instead.\n",
         "See help(aster2-package) for differences from package \"aster\"\n",
         "    and examples.\n")
}

