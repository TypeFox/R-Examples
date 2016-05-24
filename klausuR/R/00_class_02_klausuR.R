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


#' Class klausuR
#'
#' This class is used for objects that are returned by \code{\link[klausuR:klausur]{klausur}}.
#'
#' @slot results A data.frame with global results
#' @slot answ A data.frame with all given answers
#' @slot corr A vector with the correct answers
#' @slot wght A vector with the weights of items
#' @slot points A data.frame with resulting points given for the answers
#' @slot marks A vector with assignments of marks to achieved score
#' @slot marks.sum A more convenient matrix with summary information on the defined marks
#' @slot trfls A data.frame of TRUE/FALSE values, whether a subject was able to solve an item or not
#' @slot anon A data.frame for anonymous feedback
#' @slot mean A table with mean, median and quartiles of the test results
#' @slot sd Standard deviation of the test results
#' @slot cronbach Internal consistency, a list of three elements "alpha", "ci" (confidence interval 95\%) and "deleted" (alpha if item was removed)
#' @slot item.analysis A data.frame with information on difficulty, discriminant power, discriminant factor and selection index of all items.
#' @slot distractor.analysis A list with information on the selected answer alternatives for each individual item.
#' @slot test Currently an empty placeholder. Planned to hold the actual test items in future releases.
#' @slot misc Anything that was stored in the \code{misc} slot of the input data.
#' @name klausuR,-class
#' @aliases klausuR-class klausuR,-class
#' @import methods
#' @include 00_class_01_klausuR.test.R
#' @keywords classes
# @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @export
#' @rdname klausuR-class

setClass("klausuR",
    representation=representation(results="data.frame",
    answ="data.frame",
    corr="vector",
    wght="vector",
    points="data.frame",
    marks="vector",
    marks.sum="matrix",
    trfls="data.frame",
    anon="data.frame",
    mean="table",
    sd="numeric",
    cronbach="list",
    item.analysis="data.frame",
    distractor.analysis="list",
    test="klausuR.test",
    misc="data.frame"),
  prototype(results=data.frame(NULL),
    answ=data.frame(NULL),
    corr=NULL,
    wght=NULL,
    points=data.frame(NULL),
    marks=NULL,
    marks.sum=NULL,
    trfls=data.frame(NULL),
    anon=data.frame(NULL),
    mean=table(NULL),
    sd=numeric(),
    cronbach=list(alpha=NULL, ci=NULL, deleted=NULL),
    item.analysis=data.frame(NULL),
    distractor.analysis=list(NULL),
    test=new("klausuR.test"),
    misc=data.frame(NULL))
)

#setValidity("klausuR", function(object){
#  
#})
