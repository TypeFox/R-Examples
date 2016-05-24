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

#' Summary of a \code{lmmsde} Object
#' 
#' Summarizes the \code{lmmsde} object returned by the \code{\link{lmmsDE}} method. Including the models fitted, parameter used and the number of features declared as differentially expressed.
#' @param object An object of class \code{lmmsde} .
#' @param ... Additional arguments which are passed to \code{summary}.
#' @return summary of the \code{lmmsde} object.
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup) 
#' lmmsDEtest <-lmmsDE(data=kidneySimTimeGroup$data,time=kidneySimTimeGroup$time,
#'              sampleID=kidneySimTimeGroup$sampleID,group=kidneySimTimeGroup$group)
#' summary(lmmsDEtest)}

#' @method summary lmmsde
#' @export
summary.lmmsde <-function(object, ...){
            cat('Differential Expression using Linear Mixed-Effect Model Splines. \n ')
            di <- dim(object@DE)
            cat(paste('Differential expression of type "',object@type,'" were modelled for',di[1],'features.\n'))
            cat(' \n ')
            cat('Basis: \n ')
            print(object@basis)
            cat(' \n ')
            cat('Knots: \n ')
            print(object@knots)
            cat(' \n ')
            cat('Table of number of differential expressed molecules after correction for multiple testing and adj.p <=0.05.\n')
            if(di[2]>3){
              print( apply(object@DE[,seq(3,di[2],by=2)],2,function(x)sum(na.omit(x<0.05))))
              cat(' \n ')
              cat('Table of models used to model profiles: ')
              cat(' \n ')
              cat('Time: \n ')
              print(table(object@modelsUsed[,1]))
              cat(' \n ')
              cat('Group: \n ')
              print(table(object@modelsUsed[,2]))
              cat(' \n ')
              cat('Group*Time:')
              print(table(object@modelsUsed[,3]))

            }else{
              print(sum(na.omit(object@DE[,di]<0.05)))
              cat('Table of models used to model profiles: ')
              cat(' \n ')
              cat(paste(colnames(object@DE)[2],': \n '))
              print(table(object@modelsUsed))
            }

          }
