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


#' Show methods for S4 objects of classes klausuR and klausuR.mult
#'
#' Print a nice summary of the slots \code{results}, \code{anon}, \code{cronbach} and \code{item.analysis}.
#' If object is of class klausuR.mult, the global results for tests with several test forms are evaluated.
#'
#' @param object An object of class \code{klausuR} or \code{klausuR.mult}
#' @aliases show,-methods show,klausuR-method show,klausuR.mult-method
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @seealso \code{\link[klausuR:klausur]{klausur}}, \code{\link[klausuR:klausur.mufo]{klausur.mufo}}
#' @keywords methods
#' @examples
#' \dontrun{
#' klausur(data.obj)
#'
#' # multiple test forms
#' klausur.mufo(data.obj, marks=notenschluessel)
#' }
#' @export
#' @docType methods
#' @rdname show-methods
setMethod("show", signature(object="klausuR"), function(object){
  if(length(object@results) == 0){
    return(invisible(NULL))
  } else {}

  if(nrow(object@item.analysis) > 0){
    item.analysis <- data.frame(
    SD=round(object@item.analysis$Sample.SD, 2),
    Diffc=round(object@item.analysis$Difficulty, 2),
    DiscrPwr=round(object@item.analysis$Item.total, 2),
    PartWhole=round(object@item.analysis$Item.Tot.woi, 2),
    SelectIdx=round(object@item.analysis$selIdx, 2),
    Discrim=round(object@item.analysis$Discrimination, 2))
    if(length(object@item.analysis$alphaIfDeleted) > 1 && !is.na(object@item.analysis$alphaIfDeleted)){
      item.analysis["alphaIfDeleted"] <- round(object@item.analysis$alphaIfDeleted, 2)
    } else {}
    dimnames(item.analysis)[[1]] <- dimnames(object@item.analysis)[[1]]
    show.itan <- TRUE
  } else {
    show.itan <- FALSE
  }

  if(length(object@distractor.analysis) > 0){
    distractor.analysis <- object@distractor.analysis
    show.distan <- TRUE
  } else {
    show.distan <- FALSE
  }

  if(!is.null(object@cronbach$alpha)){
    cr.alpha <- paste("\t",round(object@cronbach$alpha, 2), "\n\tConfidence interval:\t",
    round(object@cronbach$ci$LCL, 2),"-",
    round(object@cronbach$ci$UCL, 2)," (95%)",sep="")
    if(!show.itan){
      cr.deleted <- round(object@cronbach$deleted, 2)
    } else {}
    show.alpha <- TRUE
  } else {
    show.alpha <- FALSE
  }

  global.results <- object@results
  global.results[["Points"]] <- round(global.results[["Points"]], digits=2)
  anon.results <- object@anon
  anon.results[["Points"]] <- round(anon.results[["Points"]], digits=2)

  cat("\nKlausuR results:")
  cat("\n\nMarks defined:\n")
  print(object@marks.sum)
  cat("\n\nGlobal results:\n")
  print(global.results)
  cat("\n\nAnonymised results:\n")
  print(anon.results)
  cat("\n\nDescriptive statistics:\n")
  print(object@mean)
  cat("\n  Sd:",round(object@sd, 2))
  if(show.alpha){
    cat("\n\nInternal consistency:\n")
    cat("\tCronbach's alpha:", cr.alpha)
    if(show.alpha && !show.itan){
      cat("\n\n")
      print(cr.deleted)
    } else {}
  } else {}
  if(show.itan){
    cat("\n\nItem analysis:\n")
    print(item.analysis)
  } else {}
  if(show.distan){
    cat("\n\nDistractor analysis:\n")
    print(distractor.analysis)
  } else {}
  cat("\n\n")
})

#' @rdname show-methods
setMethod("show", signature(object="klausuR.mult"), function(object){
  if(is.null(object@forms)){
    return(invisible(NULL))
  } else {}

  show(object@results.glob)
  cat("\n\nThese are combined results for the test forms:\n")
  print(object@forms)
  cat("\n\nYou can get the partial results for each test form by:\n")
  cat("<object.name>@results.part\n\n")
})
