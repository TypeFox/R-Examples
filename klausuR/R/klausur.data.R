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


#' A function to create data objects with given and correct answers to a test.
#'
#' \code{klausur.data} automatically parses the variable names in \code{answ}to decide \strong{which variables are actual test items},
#' if they are named according to the given scheme \code{Item###}. To help in constructing a data.frame with correct column names one can call the
#' \code{\link[klausuR:klausur.gen]{klausur.gen}} utility to generate an empty data object of a given number of items and test subjects.
#'
#' If you have \strong{items with multiple correct answers} you can easily code these as one single item: All alternatives a subject has marked should be combined
#' to a single value without spaces. The vector with correct answers will have to be coded accordingly, of course. An example: If someone marked the first,
#' third and fourth answer, you would code this as "134". See \code{\link[klausuR:klausur.gen.corr]{klausur.gen.corr}} for a helpful function to create such an
#' answer vector. Internally \code{klausur} checks for equality of given answers and correct values, that is,
#' it will only give that person a point if the correct answer was coded as "134" as well.
#'
#' \strong{Data for (Number Right) Elimination Testing}
#'
#' If your test is to be evaluated according to elimination testing (ET), number right elimination testing (NRET) or number right (NR, which is actually
#' multiple choice) scoring, the data has to be in a different format: In contrast to the usual MC procedure, ET items are answered
#' by eliminating all alternatives a subject considers \emph{wrong}; in an NRET test subjects are asked to eliminate all wrong alternatives
#' \emph{and} mark the one they consider the correct answer. That is, for both scoring functions, you need to know for each answer alternative whether
#' a subject saw it as right, wrong or was not sure and left it open.
#'
#' In this implementation, these answers are to be coded as a plus sign "\code{+}" (right answer), a minus sign "\code{-}" (wrong answer) or a zero
#' "\code{0}" (missing). If you need to code errors (like both "right" and "wrong" have been marked),use the asterisk "\code{*}" for these cases.
#' All answers to \strong{one item} belong into \strong{one column}. E.g., if you have four answer alternatives, a subject thought the second one to be the correct
#' answer and eliminated the rest, you'd have to code this item as "\code{-+--}". The same is true for the vector of correct answers, of course.
#'
#' \strong{Marks}
#'
#' The \strong{assigned marks} are expected to be in a certain format as well, as long as you don't want \code{klausur} to suggest them itself.
#' Just create an empty vector to start with (say \code{your.marks <- c()}) and fill it according to the scheme \code{your.marks[<points from>:<points to>] <- <mark>}.
#' For example: Should one get a 1.7 if in sum 27 to 30 points were achieved, you'd assign these points as indices to the vector with
#' \code{your.marks[27:30] <- "1.7"} (see example section below). It is crucial to assign marks to the whole range of points that can be achieved in the test.
#' On the other hand, it's irrelevant wheter you assign decimal marks as in the example, only integer values, a 15 marks scheme or whatever. The convenience
#' function \code{\link[klausuR:klausur.gen.marks]{klausur.gen.marks}} can assist you in creating such a valid vector.
#'
#' @param answ A \code{\link{data.frame}} which has to include at least these variables:
#'  \code{No}, \code{Name}, \code{FirstName}, \code{MatrNo}, as well as \code{Pseudonym} (optional)
#'  and variables for the answered items (according to the scheme \code{Item###},
#'  where ### is a number with leading zeros, if needed).
#' @param corr A vector with the correct answers to all items in \code{answ} (named also according to \code{Item###}).
#' @param items Indices of a subset of variables in \code{answ} to be taken as items.
#' @param marks A vector assigning marks to points achieved (see details). Leave \code{NULL} if not available.
#' @param wght A vector with weights for each item (named also according to \code{Item###}). Leave \code{NULL} if not available.
#' @param corr.key If test has several test forms: A data.frame or matrix indicating the positions of all items (columns) in all
#'    forms (rows). Must have a column called \code{Form} (like \code{answ}), and the item columns must follow the explained name
#'    scheme \code{Item###}. \code{NULL} if not needed.
#' @param rename A named vector defining if variables in \code{answ} need to be renamed into the klausuR name scheme. Accepts elements
#'    named \code{No}, \code{Name}, \code{FirstName}, \code{MatrNo}, \code{Pseudonym} and \code{Form}. The values of these elements
#'    represent the variable names of the input data.
#' @param dummies A vector of dummy variables to be created, e.g. if you don't need/want actual data in the \code{id} slot.
#'    Can include \code{"No"}, \code{"Name"}, \code{"FirstName"}, \code{"MatrNo"} and \code{"Pseudonym"}. Columns will just be filled
#'    with increasing integers.
#' @param disc.misc Logical. If \code{TRUE}, left over columns from \code{answ} will not be stored in slot \code{misc} but silently discarded.
#' @param na.rm Logical, whether cases with NAs should be ignored in \code{answ}. Defaults to TRUE.
#' @param item.prefix A named character vector with two optional elements, \code{item} and \code{corr}, defining the name prefix
#'    used for the items in the test data and the vector with correct answers, respectively. Defaults to \code{item="Item"} and \code{corr="Item"}.
#' @param sort.by A character string naming the variable to sort the \code{answ} data by. Set to \code{c()} to skip any re-ordering.
#' @param maxp Optional numeric value, if set will be forced as the maximum number of points achievable. This should actually not be needed,
#'    if your test has no strange errors. But if for example it later turns out you need to adjust one item because it has two instead of
#'    one correct answers, this option can become handy in combination with "partial" scoring and item weights.
#' @param wrong If you want full pick-n scoring: A vector similar to \code{corr}, but this time listing all alternatives that are wrong.
#' @param keep.cases A vector of \code{MatrNo} values, if you want to prevent these cases from being dropped even if they contain missing data.
#'    If not \code{NULL}, missing values in all test items are replaced by the value given to \code{recode.na}, before \code{na.rm} is evaluated.
#' @param recode.na A value to replace missing data with in all cases specified by \code{keep.cases}. Ignored if \code{keep.cases=NULL}.
#' @return An object of class \code{\link[klausuR]{klausuR.answ-class}}.
#' @export
#' @examples
#' data(antworten)
#'
#' # vector with correct answers:
#' richtig <- c(Item01=3, Item02=2, Item03=2, Item04=2, Item05=4,
#'  Item06=3, Item07=4, Item08=1, Item09=2, Item10=2, Item11=4,
#'  Item12=4, Item13=2, Item14=3, Item15=2, Item16=3, Item17=4,
#'  Item18=4, Item19=3, Item20=5, Item21=3, Item22=3, Item23=1,
#'  Item24=3, Item25=1, Item26=3, Item27=5, Item28=3, Item29=4,
#'  Item30=4, Item31=13, Item32=234)
#'
#' # vector with assignement of marks:
#' notenschluessel <- c()
#' # scheme of assignments: marks[points_from:to] <- mark
#' notenschluessel[0:12]  <- 5.0
#' notenschluessel[13:15] <- 4.0
#' notenschluessel[16:18] <- 3.7
#' notenschluessel[19:20] <- 3.3
#' notenschluessel[21]    <- 3.0
#' notenschluessel[22]    <- 2.7
#' notenschluessel[23]    <- 2.3
#' notenschluessel[24]    <- 2.0
#' notenschluessel[25:26] <- 1.7
#' notenschluessel[27:29] <- 1.3
#' notenschluessel[30:32] <- 1.0
#'
#' # now combine all test data into one object of class klausur.answ
#' data.obj <- klausur.data(answ=antworten, corr=richtig, marks=notenschluessel)
#'
#' # if that went well, get the test results
#' klsr.obj <- klausur(data.obj)

klausur.data <- function(answ, corr, items=NULL, marks=NULL, wght=NULL, corr.key=NULL, rename=c(), dummies=c(),
  disc.misc=FALSE, na.rm=TRUE, item.prefix=c(), sort.by="Name", maxp=NULL, wrong=NULL, keep.cases=NULL, recode.na=0){

  # check for var names to use
  item.prefix <- check.prefixes(prefixes=item.prefix, package="klausuR")

  # in case no items were specified, take variables of names "Item##" as items
  if(is.null(items)){
    items <- grep(paste("^(", item.prefix[["item"]], ")([[:digit:]]{1,3})$", sep=""), names(answ), ignore.case=TRUE)
  } else{}
  # we just allowed lowercase names, force these into expected diction
  names(answ)[items] <- gsub("item", item.prefix[["item"]], names(answ)[items], ignore.case=TRUE)
  names(corr) <- gsub("item", item.prefix[["corr"]], names(corr), ignore.case=TRUE)

  vars.to.rename <- names(rename)
  id.names <- c("No", "Name", "FirstName", "MatrNo")
  id.possible.names <- c(id.names, "Pseudonym", "Form")

  # check for name collisions -- do we end up with doubled varnames?
  if(length(vars.to.rename) > 0){
    if(any(vars.to.rename %in% names(answ))){
      double.vars <- vars.to.rename[vars.to.rename %in% names(answ)]
      warning(paste("Probably duplicate variable names found, please double check the outcome:\n ", paste(double.vars, collapse=", ")), call.=FALSE)
    } else {}
  } else {}
  id.invalid.names <- vars.to.rename[!vars.to.rename %in% id.possible.names]
  if(length(id.invalid.names) > 0){
    stop(simpleError(paste("Invalid variable names in 'rename':\n ",
      paste(id.invalid.names, collapse=", "))))
  } else {}

  # rename columns, if any
  for (ren.var in vars.to.rename){
    ren.from <- rename[ren.var]
    ren.to  <- ren.var
    dimnames(answ)[[2]][dimnames(answ)[[2]] == rename[ren.var]] <- ren.var
  }

  # create dummy values, if demanded
  invalid.dummies <- dummies[!dummies %in% c(id.names, "Pseudonym")]
  if(length(invalid.dummies) > 0){
    stop(simpleError(paste("Invalid variable names in 'dummies':\n ",
      paste(invalid.dummies, collapse=", "))))
  } else {
    # create dummies, if any
    for(dummy in dummies){
      answ[[dummy]] <- 1:dim(answ)[[1]]
    }
  }

  sane.data <- data.check.klausur(answ=answ, corr=corr, items=items, na.rm=na.rm, prefixes=item.prefix,
    keep.cases=keep.cases, recode.na=recode.na)
  stopifnot(scoring.check.klausur(corr=corr, marks=marks, wght=wght, score="solved", maxp=maxp))
  answ <- sane.data$answ
  items <- sane.data$items

  # convert probable factors to character, and trimming values
  found.vars <- names(answ)[names(answ) %in% c("Name", "FirstName", "Pseudonym")]
  for (char.var in found.vars){
    answ[[char.var]] <- gsub("(^[[:space:]]+)|([[:space:]]+$)", "", as.character(answ[[char.var]]))
  }
  found.vars <- names(answ)[names(answ) %in% c("No", "MatrNo")]
  for (char.var in found.vars){
    answ[[char.var]] <- as.numeric(as.character(answ[[char.var]]))
  }

  # re-order cases?
  if(length(sort.by) > 0){
    if(!sort.by %in% names(answ)){
      stop(simpleError(paste("Can't sort by '",sort.by,"', there's no such variable!", sep="")))
    } else {}
    new.order <- order(answ[[sort.by]])
    answ <- answ[new.order,]
  } else {}

  # prepare columns for resulting data.frames
  if("Pseudonym" %in% names(answ)){
    id.pseudonym <- answ[["Pseudonym"]]
  } else {
    id.pseudonym <- NA
  }
  if("Form" %in% names(answ)){
    id.form <- factor(answ[["Form"]])
  } else {
    id.form <- NA
  }
  if(isTRUE(disc.misc)){
    misc.data <- data.frame(MatrNo=answ[["MatrNo"]])
  } else {
    # collect the rest for the 'misc' slot
    unused.stuff <- answ[, !names(answ) %in% c(id.possible.names, names(answ[, items]))]
    misc.data <- data.frame(MatrNo=answ[["MatrNo"]], unused.stuff)
  }

  # force items into klausuR name scheme
  if(!identical(item.prefix[["item"]], "Item")){
    old.digits <- as.numeric(gsub(paste("^(", item.prefix[["item"]], ")([[:digit:]]{1,3})$", sep=""), "\\2", items, perl=TRUE))
    new.item.names <- gen.item.names(old.digits, prefix="Item")
    dimnames(answ[, items]) <- list(NULL,new.item.names)
  } else {}

  # create resulting object
  results <- new("klausuR.answ",
    corr=list(corr=corr, corr.key=corr.key, wrong=wrong),
    id=data.frame(
      No=answ[["No"]],
      Name=answ[["Name"]],
      FirstName=answ[["FirstName"]],
      MatrNo=answ[["MatrNo"]],
      Pseudonym=id.pseudonym,
      Form=id.form, stringsAsFactors=FALSE),
    items=data.frame(MatrNo=answ[["MatrNo"]], answ[, items], stringsAsFactors=FALSE),
    score=list(marks=marks, wght=wght, maxp=maxp),
    misc=misc.data)

  return(results)
}
