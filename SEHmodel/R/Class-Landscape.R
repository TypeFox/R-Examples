# Class-Landscape.R 
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

#' @title Class Landscape
#'
#' @description Class \code{Landscape} defines a landscape
#'
#' @details Landscape objects can be created by calling of the allocator new("Landscape", ...), or (preferred) by calling to the function \code{simulateLandscape} or \code{loadLandscape}.
#'
#' @seealso \code{simulateLandscape} , \code{loadLandscape} , \code{loadLandscapeSIG}
#'   
#' @name Landscape-class
#' @rdname Landscape-class
#' @slot thelandscape A SpatialPolygonsDataFrame
#' @slot xmin Left x-axis coordinate (in meters)
#' @slot xmax Right x-axis coordinate (in meters)
#' @slot ymin Bottom y-axis coordinate (in meters)
#' @slot ymax Top y-axis coordinate (in meters)
#' @slot n Number of fields in the landscape
#'
#' @import methods
#' @import stats
#' @import rgeos
#' @import sp
#' @importFrom grDevices colorRampPalette heat.colors rgb
#' @importFrom graphics legend lines mtext par points rect title
#' @importFrom fields Exponential
#' @importFrom MASS mvrnorm
#' @importFrom deldir deldir
#' @importFrom deldir tile.list
#' @importFrom rgdal readOGR
#'  
#' @exportClass Landscape
setClass(
  Class="Landscape",
  slots=c(thelandscape="SpatialPolygonsDataFrame",
          xmin="numeric",
          xmax="numeric",
          ymin="numeric",
          ymax="numeric",
          n="numeric"),
  prototype=list(thelandscape=new(Class="SpatialPolygonsDataFrame"),
                 xmin=0,
                 xmax=5000,
                 ymin=0,
                 ymax=5000,
                 n=1
                ),
  validity=function(object) {
    if( xmin >= xmax | ymin >= ymax) {
      stop("ERROR : [Landscape:validity] landscape coordinates error")
    }
    return(TRUE)
  }
  )

