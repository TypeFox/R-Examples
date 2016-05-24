# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Moleculesral Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Moleculesral Public License for more details.
#
# You should have received a copy of the GNU Moleculesral Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#' Summary of a \code{lmmspline} Object
#' 
#' Summarizes the \code{lmmspline} object returned by the \code{\link{lmmSpline}} method. Including the models fitted and parameter used.
#' @param object An object of class \code{lmmspline}.
#' @param ... Additional arguments which are passed to \code{summary}.
#' @return Summary of the \code{lmmspline} object.
#' @examples
#' \dontrun{
#' data(kidneySimTimeGroup)
#' # running for samples in group 1
#' G1 <- which(kidneySimTimeGroup$group=="G1")
#' testLMMSplineTG<- lmmSpline(data=kidneySimTimeGroup$data[G1,],
#'                   time=kidneySimTimeGroup$time[G1],
#'                   sampleID=kidneySimTimeGroup$sampleID[G1])
#' summary(testLMMSplineTG)}

#' @method summary lmmspline
#' @export
summary.lmmspline <- function(object,...){
            cat('Data-driven Linear Mixed-Effect Model Splines \n ')
            di <- dim(object@predSpline)
            cat(paste('Profiles were modelled for',di[1],'features with ',di[2],'time points.\n'))
            cat(' \n ')
            cat('Basis: \n ')
            print(object@basis)
            cat(' \n ')
            cat('Knots: \n ')
            cat(' \n ')
            print(object@knots)
            cat(' \n ')
            cat('Time points: \n ')
            cat(' \n ')
            print(as.numeric(colnames(object@predSpline)))
            cat(' \n ')
            cat('Table of models used to model profiles: ')
            print(table(object@modelsUsed))
            cat(' \n ')
            cat('Profiles not modelled: \n')
            print(object@errorMolecules)
          }


