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


#' A function to generate a vector with correct answers
#' 
#' Create a vector of correct answers to be used by \code{\link[klausuR:klausur]{klausur}}.
#'
#' By default answers are expected to be numeric values. You can change that to character with \code{items.char=TRUE}.
#'
#' The parameter \code{answ} is quite versatile. You can just feed it your observation data, if it complies
#' with the naming scheme for items (\code{Item###}, see also \code{\link[klausuR:klausur.gen]{klausur.gen}}), and
#' \code{klausur.gen.corr} will use all of its items automatically. Or you assign the number of items directly
#' as an integer value. If you leave \code{answ=NULL}, you will be asked for the number of items.
#' 
#' @param answ Either an object with item names in klausuR scheme (see \code{\link[klausuR:klausur]{klausur}}),
#'        e.g. your observation data, or an integer representing the maximum score of the test. If NULL, you will
#'        be asked for the maximum score.
#' @param items.char Logical, will the answers be coded as characters or integer numbers (default)?
#' @param test.forms An integer value specifying how many parallel test forms are available
#' @return A numeric or character vector (depending on the parameter \code{items.char}).
#' @keywords utilities
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @examples
#' \dontrun{
#' richtig <- klausur.gen.corr(answ=antworten)
#' }
#' @seealso \code{\link[klausuR:klausur]{klausur}}
#' @export

klausur.gen.corr <- function(answ=NULL, items.char=FALSE, test.forms=1){

    if(!is.null(answ)){
      if(is.double(answ) && length(answ) == 1)
  items <- gen.item.names(answ)
      else
  items <- names(answ[grep("Item([[:digit:]]{1,3})",names(answ))])
    }
    else {
      items <- gen.item.names(as.numeric(readline(paste("How many items are in your test?"))))
    }
  # stop if still not useful
  if(!is.vector(items))
    stop(simpleError(paste("Illegal value for items:\n",items,"\n")))
  else
    cat("Items used:\n", items,"\n")

  ## function read.answers
  # this function is called to ask for the answer to each item
  read.answers <- function(itemname, correct.answ){
    answ.item <- ""
    while(identical(answ.item, "")){
      print(correct.answ)
      answ.item <- readline(paste("Please input the correct answer for ",as.character(itemname),":", sep=""))
      if(!identical(answ.item, "")){
    if(items.char)
      correct.answ[itemname] <- as.character(answ.item)
    else
      correct.answ[itemname] <- as.numeric(answ.item)
    if(is.na(correct.answ[itemname]))
      stop(simpleError(paste("Illegal value for ",itemname,": \"",answ.item,
      "\"\nklausur.gen.corr() expects numeric values by default.\nYou probably need to use the parameter \"items.char=TRUE\"\n", sep="")))
    else
      cat("OK: Assigned answer \"",correct.answ[itemname],"\" to ",itemname,"\n", sep="")
      print(correct.answ)
      } else{}
    }
  return(correct.answ)
  } ## end function read.answers

  ## function read.test.form
  # this function is called to ask for the corrsponding item numbers or names of parallel test forms
  read.test.form <- function(itemname, form, correct.answ){
    answ.item <- ""
    while(identical(answ.item, "")){
      answ.item <- readline(paste("Form,",form,", ",as.character(itemname)," -- Please input the corresponding item in Form 1 (name or number):", sep=""))
      if(!identical(answ.item, "")){
    if(items.char)
      correct.answ[correct.answ[correct.answ$Form == form],itemname] <- as.character(answ.item)
    else
      correct.answ[correct.answ[correct.answ$Form == form],itemname] <- as.numeric(answ.item)
    if(is.na(itemname))
      stop(simpleError(paste("Illegal value for ",itemname,": \"",answ.item,
      "\"\nklausur.gen.corr() expects numeric values by default.\nYou probably need to use the parameter \"items.char=TRUE\"\n", sep="")))
    else
      cat("OK: Assigned answer \"",correct.answ[itemname],"\" to ",itemname,"\n", sep="")
      } else{}
    }
  return(correct.answ)
  } ## end function read.test.forms

  # check number of forms and create either a vector or a data.frame
  if(is.double(test.forms) && length(test.forms) == 1){
    # create empty object for the correct answers of the first test form
    correct.answ <- c()
    correct.answ.form1 <- mapply(read.answers, items, MoreArgs=list(correct.answ=correct.answ))
    correct.answ["Form"] <- factor(1, levels=1:test.forms)
    correct.answ <- rbind(correct.answ, correct.answ.form1)
    # now, read in every value!

    if(test.forms > 1){
      # only if more than one test form is needed, create a data.frame
      correct.answ <- lapply(2:test.forms, function(form){rbind(data.frame(correct.answ, stringsAsFactors=FALSE),c(Form=form,rep(NA,test.forms-1)))})
      correct.answ <- lapply(2:test.forms, function(form){mapply(read.test.form, items, MoreArgs=list(form=form, correct.answ=correct.answ))})
    } else{}
  }
  else{stop(simpleError("\"test.forms\" must be an integer value!"))}

  return(correct.answ)
}
