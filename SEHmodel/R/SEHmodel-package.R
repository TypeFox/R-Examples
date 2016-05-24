# SEHmodel-package.R 
# Part of the SEHmodel package.
#
# Copyright (C) 2015        Melen Leclerc <melen.leclerc@rennes.inra.fr>
#                           Jean-Francois Rey <jean-francois.rey@paca.inra.fr>
#                           Samuel Soubeyrand <Samuel.Soubeyrand@avignon.inra.fr>
#                           Emily Walker <emily.walker@avignon.inra.fr>
#                           INRA - BioSP Site Agroparc - 84914 Avignon Cedex 9
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#' @encoding UTF-8
#' @title Spatial Exposure-Hazard Model for Exposure and Impact Assessment on Exposed Individuals
#' @description A model coupling polygon and point processes for assessing risk due to contaminant sources and their impact on exposed individuals
#' @aliases SEHmodel-package SEHmodel
#' 
#' @author Emily Walker \email{emily.walker@@avignon.inra.fr}
#' @author Jean-Francois Rey \email{jean-francois.rey@@paca.inra.fr}
#' @author Melen Leclerc \email{melen.leclerc@@rennes.inra.fr}
#' @author Samuel Soubeyrand \email{Samuel.Soubeyrand@@avignon.inra.fr}
#' @author Marc Bourotte \email{marc.bourotte@@avignon.inra.fr}
#' 
#' Maintainer: Jean-Francois REY \email{jean-francois.rey@@paca.inra.fr}
#' @docType package
#' @name SEHmodel-package
#' @details \tabular{ll}{
#'          Package: \tab SEHmodel\cr
#'          Type: \tab Package\cr
#'          Version: \tab 0.0.11\cr
#'          Date: \tab 2016-02-04\cr
#'          License: \tab GPL (>=2)\cr
#'          }
#'
#' The SEHmodel package contains functions and methods for quantifying spatio-temporal variation in contamination risk
#' around known polygon sources of contaminants, and quantifies the impact of the contaminants on the surrounding population of individuals
#' which are located in habitat areas and are susceptible to the contaminants.
#' 
#' The package implements an exposure-hazard model based on (i) tools of stochastic geometry (marked polygon and point processes)
#' for structuring the landscape and describing the location of exposed individuals,
#' (ii) a method based on a dispersal kernel describing the dissemination of contaminant particles from polygon sources,
#' and (iii) an ecotoxicological equation describing how contaminants affect individuals of the exposed population.
#'          
#' @references Emily Walker and al. , SEHmodel: Spatially-explicit model for exposure and impact assessment. Submitted.
#' @keywords model spatial survival
#' @seealso \code{\link{demo.pollen.run}}
#' @examples
#' ## Run a simulation
#' \dontrun{
#' library("SEHmodel")
#' demo.pollen.run() 
#' }
#' @import methods
NULL
