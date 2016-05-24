# Copyright 2009-2014 Meik Michalke <meik.michalke@hhu.de>
#
# This file is part of the R package klausuR.
#
# klausuR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# klausuR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with klausuR.  If not, see <http://www.gnu.org/licenses/>.


#' Class klausuR.mult
#'
#' This class is used for objects that are returned by \code{\link[klausuR:klausur.mufo]{klausur.mufo}}.
#'
#' @slot forms A vector with the names of all test forms.
#' @slot results.part A list with the partial results of each test form
#' @slot results.glob An object of class klausuR-class with overall results
#' @name klausuR.mult-class
#' @aliases klausuR.mult-class klausuR.mult,-class
#' @import methods
#' @keywords classes
# @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @export
#' @rdname klausuR.mult-class

setClass("klausuR.mult",
  representation=representation(
      forms="vector",
      results.part="list",
      results.glob="klausuR"
  ),
  prototype(
      forms=NULL,
      results.part=list(),
      results.glob=new("klausuR")
  )
)

#setValidity("klausuR.mult", function(object){
#  
#})
