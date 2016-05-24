# pollen-data.R 
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


#' @title Pollen sources emission simulation
#' @name create.pollen.sources
#' @description Simulate pollen sources emission for maize crop. The proportion of silks 
#' emitting pollen per day was given by Frédérique Angevin (Angevin et al. 2008). 
#' @param nbOfSource Number of source fields
#' @param numberOfDay Number of days of possible emission
#' @param density Plant density plant/m\eqn{^2}
#' @param pollen Pollen production grains/plant
#' @return A matrix indexed by sources ID (in rows) and by time ( in columns) whose rows give the values of pollen emission for every source.
#' 
#' @export
create.pollen.sources <- function(nbOfSource=200,numberOfDay=60,density=runif(1,7,11),pollen=rgamma(1,shape=1.6,scale=1/(2*10^-7))) {

  # courbe d'émission de pollen par jour à partir du début de floraison : données fournies par Frédérique Angevin Proportion of emitted pollen/silks, 12 jours d'émission de pollen
  prop.pollen=c(0.0165,0.066,0.1545,0.1885,0.1735,0.156,0.1159,0.067,0.0377,0.0167,0.0055,0.0022)
  
  debut.emission=sample(1:(numberOfDay-length(prop.pollen)),nbOfSource,replace=T)
  pollen.emis=matrix(rep(0,numberOfDay*nbOfSource),ncol=numberOfDay,nrow=nbOfSource)
  
  for (i in 1:nbOfSource) {
    k=1                     
    debut=debut.emission[i]   
    for (t in debut:(debut+length(prop.pollen)-1)) {
      pollen.emis[i,t]=pollen*density*prop.pollen[k]  
      k=k+1
    }
  }
  
  return(pollen.emis)
}