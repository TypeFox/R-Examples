#########################################################################/**
# @RdocPackage R.oo
#
# \description{
#  @eval "getDescription(R.oo)"
#
#  Please note that the Rdoc syntax/grammar used to convert Rdoc comments
#  in code into Rd files is not strictly defined and is modified by the
#  need of the author. Ideally, there will be a well defined Rdoc language
#  one day.
# }
#
# \section{Installation and updates}{
#   To install this package do\cr
#
#   \code{install.packages("R.oo")}
# }
#
# \section{Dependancies and other requirements}{
#   This package requires a standard \R installation and
#   the \pkg{R.methodsS3} package.
# }
#
# \section{To get started}{
#   To get started,It is very useful to understand that:
#   \enumerate{
#   \item The @see "R.methodsS3::setMethodS3"() function, which is
#     defined in the \pkg{R.methodsS3} package (used to be part of
#     \pkg{R.oo}), is nothing but a conveniency wrapper for setting
#     up S3 methods, and automatically create an S3 generic
#     function, if missing.  For more information, see the help of
#     \pkg{R.methodsS3}.
#   \item The @see "Object" class is a top-level "root" class
#     that provides support for \emph{reference variables}.
#     Any class inheriting from this class supports
#     reference variables.
#   \item The @see "Object" class is basically a wrapper around an
#     @environment, which some additional accessors etc.  It is the
#     environment data type that provides the "emulation" of
#     reference variables - the Object class structure makes
#     it easier to extends this class and adds some level of coding
#     protection.  The Object class features is not intended for
#     referencing individual elements of basic \R data types,
#     but rather for the whole variable of such.
#     For instance, you can reassign a whole matrix \code{X} part of
#     the object this way, but you cannot reassign \code{X[1,1]}
#     without creating a completely new copy.
#   }
# }
#
# \section{Further readings}{
#   For a detailed introduction to the package see [1] (part of the
#   package distribution).
# }
#
# \section{How to cite this package}{
#   Whenever using this package, please cite [1] as\cr
#
#   @howtocite "R.oo"
# }
#
# @author
#
# \section{License}{
#   The releases of this package is licensed under
#   LGPL version 2.1 or newer.
# }
#
# \references{
#  [1] @include "../incl/BengtssonH_2003.bib.Rdoc" \cr
# }
#
# \seealso{
#   People interested in \pkg{R.oo} may also be interested in
#   packages \pkg{proto} and \pkg{mutatr}.
# }
#*/#########################################################################

