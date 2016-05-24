# Run all unit tests in installed directory unitTests
# 
# Author: Renaud Gaujoux
# Creation: 26 Oct 2011
###############################################################################

if( RcppOctave:::.isPlatformCompatible() )
    pkgmaker::utest('package:RcppOctave', quiet=FALSE)
