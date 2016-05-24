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

#' \code{lmms} class a S4 superclass to extend \code{lmmspline}  and \code{lmmsde} class.
#'
#' \code{lmms} class is a superclass for classes \code{lmmspline}  and \code{lmmsde}. These classes inherit common slots.
#'
#' @slot basis  An object of class \code{character} describing the basis used for modelling.
#' @slot knots An object of class \code{numeric}, describing the boundaries of the splines. If not defined or if basis='cubic' knots are automatically estimated using Ruppert 2002 or are the design points when using 'cubic'. 
#' @slot errorMolecules Vector of class \code{character}, describing the molecules that could not be modelled.
#'
#' @name lmms-class
#' @rdname lmms-class
#' @exportClass lmms


setClass('lmms',slots=c(basis="character", knots="numeric",errorMolecules="character"))

#' \code{lmmspline} class a S4 class that extends \code{lmms} class.
#'
#' \code{lmmspline} class inherits from class \code{lmms} and extends it with three further slots: \code{predSpline}, \code{modelsUsed}, \code{models}. The class is returned when applying \code{\link{lmmSpline}} method.
#'
#' @slot predSpline A \code{data.frame} returning the fitted values for the time points of interest.
#' @slot models  A \code{list} of class \code{\link{lm}} or  \code{\link{lme}} containing the models for every molecule
#' @slot modelsUsed A \code{list} of class \code{lm} or \code{lme}, containing the models used to model the particular feature of interest. 
#' @slot derivative A \code{logical} value indicating if the derivative was calculated.
#' 
#'
#' @name lmmspline-class
#' @rdname lmmspline-class
#' @exportClass lmmspline

setClass("lmmspline",slots= c(predSpline="data.frame", modelsUsed="numeric",models="list",derivative='logical'),contains='lmms')


#'  \code{lmmsde} class a S4 class that extends \code{lmms} class.
#'
#' \code{lmmsde} class inherits from class \code{lmms} and extends it with six further slots: DE, model.time, model.group, model.time.group, type and experiment. The class \code{lmmsde} is returned when applying \code{\link{lmmsDE}} method.
#'
#'
#' @slot DE A \code{data.frame} returning p-values and adjusted p-values using Benjamini-Hochberg correction for multiple testing of the differential expression testing over time, group or their interaction.
#' @slot modelsUsed A \code{list} of \code{lme}, containing the models used to model the particular condition of interest.  
#' @slot predTime A \code{matrix} returning the predicted time fit.
#' @slot predGroup A \code{matrix} returning the predicted group fit.
#' @slot predTimeGroup A \code{matrix} returning the predicted time group interaction fit.
#' @slot modelTime A \code{list} of class\code{\link{lme}}, containing the models for every molecule modelling the time effect. 
#' @slot modelGroup A \code{list} of class \code{\link{lme}}, containing the models for every molecule modelling group effect. 
#' @slot modelTimeGroup A \code{list} of class \code{\link{lme}}, containing the models for every molecule modelling time and group interaction effect. 
#' @slot type An object of class \code{character}, describing the test performed. 
#' @slot experiment An object of class \code{character} describing the model used to perform differential expression analysis.
#' 
#' @name lmmsde-class
#' @rdname lmmsde-class
#' @exportClass lmmsde
setClassUnion("matrixOrnumeric",c('matrix','numeric'))
setClass("lmmsde",slots= c(DE="data.frame", modelsUsed="matrixOrnumeric", predTime="matrix",predGroup="matrix",predTimeGroup="matrix",modelTime="list", modelGroup="list", modelTimeGroup="list",type="character",experiment="character"),contains='lmms')

#'  \code{noise} S4 class
#'
#' The class \code{noise} is returned when applying \code{\link{investNoise}} method.
#'
#' @slot name  \code{character} vector. The name of the trajectory. 
#' @slot RT A \code{numeric} vector, containing the time to molecule standard deviation ratios for every trajectory. 
#' @slot RI A \code{numeric} vector, containing the individual to molecule standard deviation ratios for every trajectory. 
#' @slot propMissing A \code{numeric} vector, containing the proportion of missing values for every trajectory. 
#' @slot foldChange A \code{numeric} vector, containing the maximum fold change of the mean between any two time points. 
#'  
#'
#' @name noise-class
#' @rdname noise-class
#' @exportClass noise

setClass("noise",slots= c(name="character",RT="numeric", RI="numeric", propMissing="numeric", foldChange="numeric"))


