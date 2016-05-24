#' Modelling Ordinary Differential Equations Leading to Luminescence
#'
#' A collection of function to simulate luminescence signals in the mineral quartz based on
#' published models.
#'
#' \tabular{ll}{ Package: \tab RLumModel\cr Type: \tab Package\cr Version:
#' \tab 0.1.1\cr Date: \tab 2016-05-02 \cr License: \tab GPL-3\cr }
#'
#' @name RLumModel-package
#' @docType package
#' @author \bold{Authors}
#'
#' \tabular{ll}{Johannes Friedrich \tab University of Bayreuth, Germany \cr
#' Sebastian Kreutzer \tab IRAMAT-CRP2A, Universite Bordeaux Montaigne, France\cr
#' Christoph Schmidt \tab University of Bayreuth, Germany
#'}
#'
#' \bold{Supervisor}
#'
#' Christoph Schmidt, University of Bayreuth, Germany\cr
#' Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne, France
#'
#' \bold{Support contact}
#'
#' \email{developers@@model.r-luminescence.de}\cr
#'
#' \bold{Project source code repository}\cr
#' \url{https://github.com/R-Lum}\cr
#'
#' \bold{Related projects}\cr
#' \url{http://www.r-luminescence.de}\cr
#' \url{http://cran.r-project.org/package=Luminescence}\cr
#' \url{http://shiny.r-luminescence.de}\cr
#' \url{http://cran.r-project.org/package=RLumShiny}\cr
#'
#' \bold{Package maintainer}
#'
#' Johannes Friedrich, University of Bayreuth, Germany
#' \cr \email{johannes.friedrich@@uni-bayreuth.de}
#'
#' \bold{Acknowledgement}
#'
#'  The work of Johannes Friedrich is gratefully supported by the DFG in framework of the project
#'  'Modelling quartz luminescence signal dynamics relevant for dating and dosimetry' (SCHM 305114-1)
#'
#' @keywords package
#'
#' @import Luminescence deSolve methods utils
#' @importFrom stats setNames
NULL


#' Example data (TL curve) simulated from Bailey (2001 ,fig. 1)
#'
#' @format A RLum.Analysis object containing one TL curve as RLum.Data.Curve.
#'
#' @references
#'
#' Bailey, R.M., 2001. Towards a general kinetic model for optically and thermally stimulated
#' luminescence of quartz. Radiation Measurements 33, 17-45.
#'
#' @source \bold{model_LuminescenceSignals()}
#'
#' @note This example has only one record (TL). The used sequence was
#' sequence <- list(IRR = c(temp = 20, dose = 10, DoseRate = 1),
#'                  TL = c(temp_begin = 20, temp_end = 400, heating_rate = 5))
#'
#' @keywords datasets
#' @docType data
#' @aliases model.output
#' @examples
#'
#' data(ExampleData.ModelOutput)
#' TL_curve <- get_RLum(model.output, recordType = "TL$", drop = FALSE)
#'
#' ##plot TL curve
#' plot_RLum(TL_curve)
#'
#' TL_concentrations <- get_RLum(model.output, recordType = "(TL)", drop = FALSE)
#' plot_RLum(TL_concentrations)
#'
#'
#' @name ExampleData.ModelOutput
NULL
