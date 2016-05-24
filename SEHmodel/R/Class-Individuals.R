# Class-Individuals.R 
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

#' @title Class Individuals
#' @name Individuals-class
#' 
#' @description
#' \code{Individuals} Class consists of spatio-temporal parameters about exposed populations.
#' 
#' Individual gets coordinates (as SpatialPoints), a date of birth, a life duration, an intern toxic concentration over the time and a toxic threshold (max value of toxic before death).
#' 
#' Each individual in an Individuals object is identified by an ID which is used as index to access attributes in the object.
#' 
#' @details Objects can be created by calling of the allocator new("Individuals", ...), or (preferred) by calling one of the wrapped functions \link{simulateIndividuals} or \link{loadIndividuals}.
#'
#' @seealso \code{\link{simulateIndividuals}} , \code{\link{loadIndividuals}}
#'
#' @rdname Individuals-class
#' 
#' @slot n individuals Number
#' @slot coordinate individuals coordinates (as SpatialPoints)
#' @slot xmin x-axis left value
#' @slot xmax x-axis right value
#' @slot ymin y-axis bottom value
#' @slot ymax y-axis top value
#' @slot dob Date of birth (as vector)
#' @slot life_duration individuals life duration (as vector)
#' @slot intern_toxic individuals intern toxic concentration in time (as matrix)
#' @slot toxic_threshold individuals max toxic concentration leading to death (as vector)
#' @slot mintime Start simulation time
#' @slot maxtime End simulation time
#'  
#' @exportClass Individuals
setClass(Class="Individuals",
         slots=c(n="numeric",
                 coordinate="SpatialPoints",
                 xmin="numeric",
                 xmax="numeric",
                 ymin="numeric",
                 ymax="numeric",
                 dob="vector",
                 #adult="vector",
                 life_duration="vector",
                 intern_toxic="matrix",
                 toxic_threshold="vector",
                 mintime="numeric",
                 maxtime="numeric"),
         prototype=list(n=0,
                        coordinate=new("SpatialPoints"),
                        xmin=0,
                        xmax=1000,
                        ymin=0,
                        ymax=1000,
                        dob=NULL,
                        life_duration=NULL,
                        intern_toxic=NULL,
                        toxic_threshold=NULL,
                        mintime=1,
                        maxtime=2),
         validity=function(object) {
           if(mintime == 0 | mintime>maxtime) {
             stop("ERROR : [Individuals:validity] mintime ",mintime," or maxtime ",maxtime," bad values")
           }
           if(n <0 & length(dob) != n & length(toxic_threshold) != n & length(intern_toxic) != n*maxtime) {
             stop("ERROR : [Individuals:validity] individuals numbers or data length not match")
           }
           if( xmin >= xmax | ymin >= ymax) {
             stop("ERROR : [Individuals:validity] landscape size error")
           }
           return(TRUE)
         }
)
