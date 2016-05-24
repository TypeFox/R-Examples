##
##  z z z . R
##


# Define 'ans' as in Matlab
# But: "Package namespaces are locked when loaded!"
# makeActiveBinding("ans", function() .Last.value, .GlobalEnv)

.pracmaEnv <- new.env()
assign("elapsedTime", 0, envir = .pracmaEnv)

.onLoad <- function(libname, pkgname) {
    # require(some_packages)

    # Load dynamic libraries
    # library.dynam(pkg, pkg, lib)

    environment(.pracmaEnv) <- asNamespace("pracma")

    # packageStartupMessage(
    #     paste("pracma Package Version 1.4.6\n",
    #           "Practical Numerical Math Functions\n",
    #           "Copyright (c) 2011-2013 Hans W Borchers\n",
    #     sep='', collapse=''))
}
