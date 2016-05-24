#' Tools for Thermoluminescences Dating
#'
#' A series of functions for thermoluminescence dating using the MAAD or the
#' SAR protocol. This package adds to the R package \link{Luminescence}.
#'
#' \tabular{ll}{  Package: \tab TLdating \cr
#'                Type: \tab Package\cr
#'                Version: \tab 0.1.1 \cr
#'                Date: \tab 2016-03-01 \cr
#'                License: \tab GPL-3
#'             }
#'
#' @name TLdating-package
#' @aliases TLdating-package TLdating
#' @docType package
#'
#' @author
#'  \bold{Authors} \cr
#'  \tabular{ll}{  David Strebler, \tab University of Cologne, Germany
#'  }
#'
#'  \bold{Beta-tester} \cr
#'  \tabular{ll}{  Anja Zander, \tab University of Cologne, Germany
#'  }
#'
#'  \bold{Supervisor} \cr
#'  \tabular{ll}{  Helmut Brückner, \tab University of Cologne, Germany \cr
#'                 Dominik Brill, \tab University of Cologne, Germany
#'  }
#'
##  \bold{Support contact} \cr
##  \email{david.strebler@uni-koeln.de} \cr
#'
##  \bold{Bug reporting} \cr
##  \email{david.strebler@uni-koeln.de} \cr
#'
##  \bold{Project website} \cr
##  ... \cr
#'
#'  \bold{Project source code repository} \cr
#'  \url{https://github.com/dstreble/TLdating} \cr
#'
#'  \bold{Related package projects}\cr
#'    \url{http://www.r-luminescence.de} \cr
#'    \url{http://cran.r-project.org/package=Luminescence}\cr
#'
#'  \bold{Package maintainer} \cr
#'  David Strebler, Geographisches Institut, Universitat zu Koeln, Cologne, Germany. \cr
#'  \email{david.strebler@uni-koeln.de}
#'
#'  \bold{Acknowledgement} \cr
#'  This project is realized in the context of the CRC 806 “Our Way to Europe” (\url{http://www.sfb806.uni-koeln.de/})
#'  which is funded by the German Research foundation (DFG). \cr
#'  The code and the structure of this package is partially based on those from the R package \link{Luminescence}
#'  version 0.5.1 developed by S. Kreutzer, M. Dietze, C. Burow, M.C. Fuchs, C. Schmidt, M. Fisher
#'  and R.K. Smedley. \cr
#'
## @references
#'
#' @keywords package
#'
#' @import Luminescence methods
#' @importFrom grDevices rainbow
#' @importFrom graphics abline arrows axis layout lines mtext par plot points segments title
#' @importFrom stats lm nls nls.control smooth.spline
#' @importFrom utils glob2rx
#' @importFrom gplots textplot

NULL


