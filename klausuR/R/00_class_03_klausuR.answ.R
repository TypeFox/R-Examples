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


#' Class klausuR.answ
#'
#' This class is used for objects needed by \code{\link[klausuR:klausur]{klausur}}.
#' They contain all relevant data to calculate test results.
#'
#' @note The slots \code{}, \code{id}, \code{items} and \code{misc}, must have the same number of rows and contain copies of the colum \code{MatrNo}
#' for identification.
#
#' @slot corr Contains three elements:
#'    \itemize{
#'      \item{\code{corr}} {The correct answers to each item.}
#'      \item{\code{corr.key}} {An optional data.frame or matrix for test with multiple test forms, indicating the positions
#'        of all items (columns) in all forms (rows). Must have a column called \code{Form} (like in \code{id}), and the
#'        item columns must follow the usual name scheme \code{Item###}. \code{NULL} if not needed.}
#'      \item{\code{wrong}} {For pick-n scoring, this is the inverse of \code{corr}, i.e., all wrong item alternatives.}
#'    }
#' @slot id Contains the columns \code{No}, \code{Name}, \code{FirstName}, \code{MatrNo}, \code{Pseudonym} and \code{Form}.
#' @slot items Contains a copy of \code{id$MatrNo} and all answers to the test items (one item per column).
#' @slot score Contains three elements:
#'    \itemize{
#'      \item{\code{marks}} {The assigned marks for achieved points (\code{NULL} if none)}
#'      \item{\code{wght}} {Weights for each item (\code{NULL} if none)}
#'      \item{\code{maxp}} {Optional, to force a certain maximum points value (\code{NULL} if none)}
#'    }
#' @slot test Currently an empty placeholder. Planned to hold the actual test items in future releases.
#' @slot misc Any additional data you'd like to be stored along with \code{id} and \code{items},
#'    e.g. table data from/for other software products. Won't be used for anything.
#' @name klausuR.answ,-class
#' @aliases klausuR.answ-class klausuR.answ,-class
#' @import methods
#' @include 00_class_01_klausuR.test.R
#' @keywords classes
# @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @export
#' @rdname klausuR.answ-class

setClass("klausuR.answ",
  representation=representation(
    corr="list",
    id="data.frame",
    items="data.frame",
    score="list",
    test="klausuR.test",
    misc="data.frame"
  ),
  prototype(
    corr=list(corr=c(Item1=NA), corr.key=NULL, wrong=NULL),
    id=data.frame(No=NA, Name=NA, FirstName=NA, MatrNo=NA, Pseudonym=NA, Form=NA),
    items=data.frame(MatrNo=NA, Item1=NA),
    score=list(marks=NULL, wght=NULL, maxp=NULL),
    test=new("klausuR.test"),
    misc=data.frame(MatrNo=NA)
  )
)

setValidity("klausuR.answ", function(object){
  obj.corr.l  <- object@corr
  obj.id     <- object@id
  obj.items   <- object@items
  obj.score  <- object@score
  obj.misc   <- object@misc

  id.names <- c("No", "Name", "FirstName", "MatrNo")
  id.possible.names <- c(id.names, "Pseudonym", "Form")
  
  # check if all slots have the same number of rows
  if(nrow(obj.id) != nrow(obj.items) | nrow(obj.id) != nrow(obj.misc)){
    stop(simpleError("All slots must have the same number of rows!"))
  } else {}
  # check if all slots have a column MatrNo (id will be checked later)
  if(!"MatrNo" %in% names(obj.items)){
    stop(simpleError(paste("Missing variable 'MatrNo' in slot 'items'!")))
  } else if(!"MatrNo" %in% names(obj.misc)){
    stop(simpleError(paste("Missing variable 'MatrNo' in slot 'misc'!")))
  } else {}
  if(length(obj.items$MatrNo) != length(unique(obj.items$MatrNo))){
    dupl.matn <- obj.items$MatrNo[duplicated(obj.items$MatrNo)]
    stop(simpleError(paste("Duplicate entries for 'MatrNo' are not allowed:\n", paste(dupl.matn, collapse=", "))))
  } else {}
  # everything in its place?
  if(any(!id.names %in% names(obj.id))){
      missing.id.vars <- id.names[!id.names %in% names(obj.id)]
      stop(simpleError(paste("Missing variables in slot 'id':\n ", paste(missing.id.vars, collapse=", "))))
  } else {}
  if(any(!names(obj.id) %in% id.possible.names)){
      invalid.id.vars <- names(obj.id)[!names(obj.id) %in% id.possible.names]
      stop(simpleError(paste("Invalid columns in slot 'id', put them in slot 'misc':\n ", paste(invalid.id.vars, collapse=", "))))
  } else {}
  # check names in slot items
  colnames.items <- names(obj.items)[names(obj.items) != "MatrNo"]
  items.valid.names <- grepl("^(item|Item)([[:digit:]]{1,3})$", colnames.items) # TRUE for valid names
  if(sum(items.valid.names) < 1){
    stop(simpleError(paste("You didn't define any items in slot 'items'!")))
  } else {}
  if(sum(!items.valid.names) > 0){
    stop(simpleError(paste("Invalid columns in slot 'items', put them in slot 'misc':\n ",
      paste(colnames.items[!items.valid.names], collapse=", "))))
  } else {}
  # check corr elements
  corr.names <- names(obj.corr.l)
  invalid.corr.names <- corr.names[!corr.names %in% c("corr","corr.key","wrong")]
  missing.corr.names <- corr.names[!c("corr","corr.key","wrong") %in% corr.names]
  if(length(invalid.corr.names) > 0){
    stop(simpleError(paste("Invalid elements in slot 'corr':\n ",
      paste(invalid.corr.names, collapse=", "))))
  } else {}
  if(length(missing.corr.names) > 0){
    stop(simpleError(paste("Missing elements in slot 'corr':\n ",
      paste(missing.corr.names, collapse=", "))))
  } else {}
  obj.corr <- obj.corr.l$corr
  # check if each item has it's correct answer
  if(!setequal(colnames.items, names(obj.corr))){
    miss.items.corr <- names(obj.corr)[!names(obj.corr) %in% colnames.items]
    miss.items.items <- colnames.items[!colnames.items %in% names(obj.corr)]
    stop(simpleError(paste("Please check:\n  ", paste(miss.items.corr, miss.items.items, collapse=", "), "\n  The number of items differs between observed and correct answers!", sep="")))
  }
  ## TODO: check corr.key
  obj.corr.key <- obj.corr.l$corr.key
  ## TODO: check wrong (e.g., for reoccurences from corr!)
  obj.corr.wrong <- obj.corr.l$wrong

  # check score elements
  score.names <- names(obj.score)
  invalid.score.names <- score.names[!score.names %in% c("marks","wght","maxp")]
  missing.score.names <- score.names[!c("marks","wght") %in% score.names]
  if(length(invalid.score.names) > 0){
    stop(simpleError(paste("Invalid elements in slot 'score':\n ",
      paste(invalid.score.names, collapse=", "))))
  } else {}
  if(length(missing.score.names) > 0){
    stop(simpleError(paste("Missing elements in slot 'score' (set them to NULL if empty):\n ",
      paste(missing.score.names, collapse=", "))))
  } else {}
  
  # check if MatrNo is identical in all slots
  if(any(!setequal(obj.id$MatrNo, obj.items$MatrNo), !setequal(obj.id$MatrNo, obj.misc$MatrNo))){
    stop(simpleError(paste("Values of column 'MatrNo' differ between slots, but must be identical!")))
  }

})
