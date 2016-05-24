## $Date$
## $Id$

# grab the namespace
.NAMESPACE <- environment()

.onLoad <- function(lib, pkg) {
    library.dynam("RDieHarder", pkg, lib )

    ## assign lists of generators and tests to hidden variables in namespace
    assign(".dieharder.generators", dieharderGenerators(), .NAMESPACE)
    assign(".dieharder.tests", dieharderTests(), .NAMESPACE)
}

