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


#' Generate empty data sets for the use with klausur
#' 
#' @param items An integer declaring the number of items to be created
#' @param obs Integer, numer ob observations
#' @param items.char Logical, will the answers be coded as characters or integer numbers (default)?
#' 
#' @return A data.frame containing the variables "No", "Name", "FirstName", "MatrNo", "Pseudonym", "Form" and the number of items as needed.
#' @seealso \code{\link[klausuR:klausur]{klausur}}, \code{\link[klausuR:compare]{compare}},
#'   \code{\link[klausuR:klausur.gen.marks]{klausur.gen.marks}}, \code{\link[klausuR:klausur.gen.corr]{klausur.gen.corr}}
#' @keywords datagen
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @examples
#' antworten2 <- klausur.gen(items=20,obs=40)
#' @export

klausur.gen <- function(items=NULL, obs=1, items.char=FALSE){
  # to avoid NOTEs from R CMD check:
  daten <- NULL

  if(!is.numeric(items) || (floor(items) != items) || !(items < 1000))
    stop(simpleError("'items' must be an integer < 1000"))
  if(!is.numeric(obs) || (floor(obs) != obs))
    stop(simpleError("'obs' must be an integer"))

  # generate item names and compute the needed number of leading zeros
  # gen.item.names() is an internal function of package klausuR
  item.names <- gen.item.names(items)

  data.filler <- if(items.char){as.character(c(NA))} else {as.integer(NA)}

  eval(parse(text=paste("daten <- data.frame(
      No=c(1:obs),
      Name=as.character(c(NA)),
      FirstName=as.character(c(NA)),
      MatrNo=as.integer(c(NA)),
      Pseudonym=as.character(c(NA)),
      Form=as.character(c(NA)),",
      paste(item.names, "=data.filler,", collapse=""),
      "stringsAsFactors=FALSE)")))

  return(daten)
}
