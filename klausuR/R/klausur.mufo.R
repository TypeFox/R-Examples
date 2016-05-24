# Copyright 2009-2015 Meik Michalke <meik.michalke@hhu.de>
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


#' Evaluate multiple choice tests with several test forms
#'
#' This function can be used to evaluate tests that have several test forms. Please be aware that its results only make sense
#' if each test form uses the same items, only in a different order.
#'
#' Firstly, \code{klausur.mufo} will compute partial results for each parallel form, and in the end combine these to global
#' results. Cronbach alpha and item analysis will be calculated for all subjects accordingly, therefore the test items of
#' all tests will be re-ordered to fit the order of the first given test form (this does not apply to the partial results).
#'
#' The parameters are mostly the same as those for \code{\link[klausuR:klausur]{klausur}}. However, in the \code{data} object the
#' slot \code{corr} must also contain \code{corr.key}, to communicate the order of items in each test form, and the slot \code{id}
#' needs one additional variable called \code{Form}.
#' 
#' An example: You have prepared a test in two different parallel forms "A" an "B", So in addition to the variables in \code{data@@id}
#' you need to create a variable called \code{Form}, to document which test subject was given which test form. Since form "B" holds the same
#' items as form "A", only in a different order, we only need to define these positions and we're done. Therefore \code{corr.key} must
#' be a matrix or data.frame, again with a column called "Form", one column for each item, and one row of data for each test form. That is, you'd need
#' one row for test form "A" and one for test form "B", giving an index for each item where it is placed in the form. For "A" this is
#' simply ascending numbers from 1 to how many questions you asked, but for row "B" each number indicates at which position an item
#' of "A" is to be found. See the example below.
#' 
#' @param data An object of class \code{\link[klausuR]{klausuR.answ-class}}.
#' @param marks A vector assigning marks to points achieved (see details). Alternatively, set it to \code{"suggest"} to let
#'    \code{\link[klausuR:klausur.gen.marks]{klausur.gen.marks}} calculate suggestions under the assumption of normal distribution.
#'    If \code{NULL}, this value must be set in the \code{data} object.
#' @param mark.labels If \code{marks="suggest"}, use these as the marks you want to give.
#' @param items Indices of a subset of variables in \code{answ} to be taken as items.
#' @param wght A vector with weights for each item (named also according to \code{Item###}). If \code{NULL}, the value from the \code{data} object
#'    will be used.
#' @param score Specify the scoring policy, must be one of \code{"solved"} (default), \code{"partial"}, \code{"liberal"},
#'    \code{"NR"}, \code{"ET"}, \code{"NRET"}, or \code{"NRET+"}.
#' @param matn A matriculation number of a subject, to receive detailed results for that subject.
#' @param na.rm Logical, whether cases with NAs should be ignored in \code{data}. Defaults to TRUE.
#' @param cronbach Logical. If TRUE, Cronbach's alpha will be calculated.
#' @param item.analysis Logical. If TRUE, some usual item statistics like difficulty and discriminatory power will be calculated.
#'  If \code{cronbach} is TRUE, too, it will include the alpha values if each item was deleted.
#' @return An object of class \code{\link[klausuR]{klausuR.mult-class}} with the following slots.
#'  \item{forms}{A character vector naming all test forms}
#'  \item{results.part}{A list of objects of class \code{klausuR}, holding all partial results}
#'  \item{results.glob}{An object of class \code{klausuR} with the global results}
#'  Not all slots are shown by default (refer to \code{\link[klausuR:show-methods]{show}}).
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @seealso \code{\link[klausuR:klausur]{klausur}}, \code{\link[klausuR:klausur.data]{klausur.data}}
#' @keywords misc
#' @import psychometric
#' @export
#' @examples
#' # this will create the data.frame "antworten.mufo"
#' # and the matrix "corr.key"
#' data(antworten.mufo)
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
#' mufo.data.obj <- klausur.data(answ=antworten.mufo, corr=richtig, marks=notenschluessel, corr.key=corr.key)
#' # expect some warnings here, because some items have no variance
#' # in their subtest results, hence item analysis fails on them
#' klsr.mufo.obj <- klausur.mufo(mufo.data.obj)

klausur.mufo <- function(data, marks=NULL, mark.labels=NULL, items=NULL, wght=NULL, score="solved", matn=NULL, na.rm=TRUE, cronbach=TRUE, item.analysis=TRUE){
  # to avoid NOTEs from R CMD check:
  MatrNo <- NULL

  if(!inherits(data, "klausuR.answ")){
    stop(simpleError("'data' must be of class 'klausuR.answ'!"))
  } else {}

  corr <- slot(data, "corr")$corr
  corr.key <- slot(data, "corr")$corr.key
  if(is.null(wght)){
    wght <- slot(data, "score")$wght
    if(is.null(wght)){
      # for the results, create a vector of 1's by number of items
      wght.results <- rep(1, length(corr))
    } else {
      wght.results <- wght
    }
  } else {
    wght.results <- wght
  }
  if(is.null(marks)){
    marks <- slot(data, "score")$marks
  } else {}

  ## first we'll check if the data is sane
  # are we given any parallel forms at all?
  if(is.null(slot(data, "id")[["Form"]]) || is.null(dim(corr.key)) || is.null(corr.key[, "Form"])){
    stop(simpleError("Both answ and corr.key must contain the variable \"Form\"!"))
  } else {}
  # ok, so how many test forms are there?
  test.forms.answ <- levels(as.factor(slot(data, "id")[["Form"]]))
  test.forms.corr <- levels(as.factor(corr.key[, "Form"]))
  # the relevant value is what's in answ. it won't hurt if corr.key has more,
  # as long as it holds all relevant data for the forms in answ
  missing.forms <- !is.element(test.forms.answ, test.forms.corr)
  # if they are not compatible, we'll cry a little. and then stop.
  if(sum(missing.forms) > 0){
    stop(simpleError(paste("Sorry, but corr.key lacks information on test form", paste(test.forms.answ(missing.forms), collapse=", "))))
  } else {}

  # if a user accidently(?) submits only one test form, it's clearly a case for klausur()
  # we'll just hand it over and return the results
  if(length(test.forms.answ) == 1){
    warning("Only one test form was supplied. Called klausur() instead.")
    klausur.mufo.results <- klausur(data=data, marks=marks,
            mark.labels=mark.labels, items=items, wght=wght, score=score, matn=matn,
            na.rm=na.rm, cronbach=cronbach, item.analysis=item.analysis)
    ## calculation would end here if only one form was submitted
  } else {
    ## separate forms and calculate partial results
    # we'll call klausur() for each test form, and later combine the results
    klausur.mufo.results <- new("klausuR.mult")
    # first save the names of all test forms to the first slot
    slot(klausur.mufo.results, "forms") <- test.forms.answ
    # create an empty object of class klausuR to store combined global results
    klausur.mufo.global <- new("klausuR", corr=corr, marks=marks, marks.sum=marks.summary(marks),
      wght=wght.results, misc=as.data.frame(data@misc, stringsAsFactors=FALSE))
    for(single.test in test.forms.answ){
      corr.inices <- as.numeric(corr.key[corr.key[, "Form"] == single.test,][!dimnames(corr.key)[[2]] == "Form"])
      corr.part <- corr[corr.inices]
      names(corr.part) <- names(corr)
      subj.in.test <- slot(data, "id")$Form == single.test
      data.part <- new("klausuR.answ",
        corr=list(corr=corr.part, corr.key=NULL, wrong=NULL),
        id=subset(slot(data, "id"), subj.in.test),
        items=subset(slot(data, "items"), subj.in.test),
        score=list(marks=marks, wght=wght),
        misc=subset(slot(data, "misc"), subj.in.test))
      klausur.part.results <- klausur(data=data.part, marks=marks,
              mark.labels=mark.labels, items=items, wght=wght, score=score, matn=matn,
              na.rm=na.rm, cronbach=cronbach, item.analysis=item.analysis)
      result.part <- list(Form=single.test, Results=klausur.part.results)
      slot(klausur.mufo.results, "results.part")[[single.test]] <- result.part
      # append partial results to global results
      slot(klausur.mufo.global, "results") <- rbind(slot(klausur.mufo.global, "results"), slot(klausur.part.results, "results"))
      slot(klausur.mufo.global, "answ")    <- rbind(slot(klausur.mufo.global, "answ"), klausur.reorderItems(slot=slot(klausur.part.results, "answ"), order=corr.inices))
      slot(klausur.mufo.global, "trfls")   <- rbind(slot(klausur.mufo.global, "trfls"), klausur.reorderItems(slot=slot(klausur.part.results, "trfls"), order=corr.inices))
      slot(klausur.mufo.global, "anon")    <- rbind(slot(klausur.mufo.global, "anon"), slot(klausur.part.results, "anon"))
    }
      ## combine global results
      # missing pieces we need to calculate
      slot(klausur.mufo.global, "mean")  <- summary(slot(klausur.mufo.global, "results")$Points)
      slot(klausur.mufo.global, "sd")  <- sd(slot(klausur.mufo.global, "results")$Points)
    ## psychometic quality of the items
    if(isTRUE(cronbach)){
      # calling an internal function which is
      # using alpha() from package "psychometric"
      cron.alpha.list <- calc.cronbach.alpha(subset(slot(klausur.mufo.global, "trfls"), select=-MatrNo))
    } else {
      cron.alpha.list <- list(alpha=NULL, ci=NULL, deleted=NULL)
    }
    if(isTRUE(item.analysis)){
      # calling another internal function which is also
      # using alpha() from package "psychometric"
      item.analyse <- calc.item.analysis(subset(slot(klausur.mufo.global, "trfls"), select=-MatrNo), cron.alpha.list)
      distractor.analysis <- distrct.analysis(
          answ=slot(klausur.mufo.global, "answ"),
          corr=corr,
          points=slot(klausur.mufo.global, "trfls"),
          results=data.frame(MatrNo=slot(klausur.mufo.global, "results")$MatrNo, Points=slot(klausur.mufo.global, "results")$Points),
          partWhole=isTRUE(score %in% c("solved", "NR"))
        )
    } else {
      item.analyse <- data.frame(NULL)
      distractor.analysis <- list(NULL)
    }

    slot(klausur.mufo.global, "cronbach") <- cron.alpha.list
    slot(klausur.mufo.global, "item.analysis") <- item.analyse
    slot(klausur.mufo.global, "distractor.analysis") <- distractor.analysis
    
    slot(klausur.mufo.results, "results.glob") <- klausur.mufo.global
  }

  return(klausur.mufo.results)
}
