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

#' Summary of a \code{noise} Object
#' 
#' Summarizes the \code{noise} object returned by the \code{\link{investNoise}} method.
#' @param object An object of class \code{noise}.
#' @param ... ignored
#' @return Summary of the \code{noise} object.
#' @examples
#' \dontrun{
#' data(kidneySimTimeGroup)
#' # running for samples in group 1
#' G1 <- which(kidneySimTimeGroup$group=="G1")
#' noiseTest<- investNoise(data=kidneySimTimeGroup$data[G1,],
#'                   time=kidneySimTimeGroup$time[G1],
#'                   sampleID=kidneySimTimeGroup$sampleID[G1])
#' summary(noiseTest)}

#' @method summary noise
#' @export
summary.noise <- function(object,...){
            cat('Time course noise level summary \n ')
            cat(' \n ')
            cat('R_T level info: \n ')
            cat(' \n ')
            print( summary(object@RT)[c(1,4,6)])
            cat(' \n ')
            cat('R_I level info: \n ')
            cat(' \n ')
            print( summary(object@RI)[c(1,4,6)])
            cat(' \n ')
            cat('Proportion of missing values: \n ')
            cat(' \n ')
            print( summary(object@propMissing)[c(1,4,6)])
            cat(' \n ')
            cat('Fold changes: \n ')
            cat(' \n ')
            print(summary(object@foldChange)[c(1,4,6)])
          }


