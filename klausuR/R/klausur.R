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


#' Evaluate multiple choice tests
#'
#' The function \code{klausur} expects an object of class \code{\link[klausuR]{klausuR.answ-class}}, containing some
#' identification data on all subjects and their answers to the test items, a vector with the correct answers, and optionally a vector with
#' marks assigned to the points achieved. It will compute global test results as well as some item analysis (including Cronbach's alpha,
#' discriminatory power and Lienert's selection index of the test items), and anonymous feedback for the test subjects.
#'
#' For details on the ecpected data structure refer to \code{\link[klausuR:klausur.data]{klausur.data}},
#'
#' \strong{Scoring functions}
#'
#' In combination with multiple (correct) answers for certain items you can specify one of six scoring policies via the \code{score}-parameter.
#' If you set it to something othe than \code{"solved"}, as the names may suggest you \strong{allow partially given answers} under the condition that the test subject didn't
#' check more alternatives than there are correct ones (that is, if you checked four alternatives where three correct ones were possible, you're out):
#' \itemize{
#'    \item{\code{"solved"}} {Multiple Choice: Check all correct alternatives. This is the default, which means that one will only get any points for an
#'      item if the answer was 100\% correct (that is, all or nothing).}
#'    \item{\code{"partial"}} {Multiple Choice: Check all correct alternatives, allow partially given answers, but none of the distractors must be checked.}
#'    \item{\code{"liberal"}} {Multiple Choice: Check all correct alternatives, allow partially given answers, even distractors can be checked.}
#'    \item{\code{"pick-n"}} {Multiple Choice: Check all correct alternatives, allow partially given answers, even distractors can be checked. Difference to \code{"liberal"} is
#'      that you will also get points for unchecked distractors.}
#'    \item{\code{"ET"}} {Elimination Testing: In contrast to the usual MC procedure, eliminate/strike all \emph{wrong} alternatives.}
#'    \item{\code{"NRET"}} {Number Right Elimination Testing: Like ET+MC, eliminate/strike all \emph{wrong} alternatives \emph{and} check the correct one.}
#'    \item{\code{"NRET+"}} {Number Right Elimination Testing, more strict: Like NRET, but if more alternatives are checked right than there are right anwers,
#'      it will automatically yield to 0 points for that item.}
#'    \item{\code{"NR"}} {Number Right: The usual MC scoring, but works with ET/NRET data. It was implemented for completeness, e.g. to compare results
#'      of different scoring techniques.}
#' }
#'
#' An example for \code{"solved"}, \code{"partial"} and \code{"liberal"}: If an item has five answer alternatives, the correct answer is "134" and a subject checked "15",
#' \code{"solved"} will give no point (because "15" is not equal to "134"), as will \code{"partial"} (because "5" is wrong), but \code{"liberal"} will give 1/3 (because "1" is correct),
#' and \code{"pick-n"} will give 2/5 (because "1" was correctly checked and "2" correctly unchecked).
#'
#' \strong{(Number Right) Elimination Testing}
#'
#' Note that \code{"ET"}, \code{"NRET"}/\code{"NRET+"} and \code{"NR"} will disable \code{wght} as of now,
#' and \strong{need the data in different format} than the other scoring functions (see \code{\link[klausuR:klausur.data]{klausur.data}} for details).
#' \code{klausur} will evaluate each answer individually and sum up the points for each item. The alternative-wise evaluations will be documented in the
#' \code{trfls} slot of the results.
#' Therefore, in these cases that matrix is not boolean, but more complex. For each item and each subject, a character string represents
#' the evaluated answer alternatives, with the following elements:
#' \itemize{
#'    \item{\code{P}} {True positive: Alternative was checked as right and is right. \emph{Points: +1 (+ constant)}}
#'    \item{\code{p}} {False positive: Alternative was checked as right but is wrong. \emph{Points: 0 (+ constant)}}
#'    \item{\code{N}} {True negative: Alternative was checked as wrong and is wrong. \emph{Points: +1 (+ constant)}}
#'    \item{\code{n}} {False negative: Alternative was checked as wrong but is right. \emph{Points: -(alternatives-1) (+ constant)}}
#'    \item{\code{0}} {Missing: Alternative wasn't checked at all. \emph{Points: 0 (+ constant)}}
#'    \item{\code{*}} {Error: Alternative was checked both wrong and right. \emph{Points NRET+: 0 (+ constant);
#'      NR scores 1 point if this was the correct alternative, ET 1 point if it hit a wrong one,
#'      and NRET sums up the points for both the positive and negative answer (all + constant)}}
#' }
#' An example: If we have an item with four alternatives, and the third one is right (i.e., "\code{--+-}"), and a test subject considered the first alternative
#' to be correct and eliminated all others (i.e., "\code{+---}"), it would be evaluated as "\code{pNnN}", that is 0+1-3+1=-1 point, not considering the constant.
#' As you can see, it would be possible to end up with a negative sum of points. If you consider how in the end a mark will be assigned to the achieved points,
#' this would be a problem, because a vactor cannot have negative indices. To circumvent this issue, klausuR automatically adds a constant to all results, so
#' that the worst possible result is not negative but 0. This constant is simply (alternatives-1), i.e. 3 for the example. In other words, if our test had 10
#' such items, the results minus 30 would be equivalent to scoring without that constant. You can use \code{\link[klausuR:nret.rescale]{nret.rescale}} to remove
#' the constant from the results afterwards.
#' 
#' \strong{Marks}
#'
#' The \strong{assigned marks} are expected to be in a certain format as well (see \code{\link[klausuR:klausur.data]{klausur.data}} for details),
#' as long as you don't want \code{klausur} to suggest them itself. If you want to let klausuR make a suggestion, set \code{marks="suggest"},
#' and \code{\link[klausuR:klausur.gen.marks]{klausur.gen.marks}} kicks in and takes either the \code{mark.labels} you have defined here or will
#' ask you step by step. See the documentation of that function for details. To see the suggested result in detail, have a look at the slot
#' \code{marks} of the returned object.
#'
#' To calculate Cronbach's alpha and item analysis methods from the package \code{\link[psychometric]{psychometric}} are used. Lienert's selction index
#' ("Selektionskennwert") aims to consider both discriminatory power (correlation of an item with the test results) and difficulty to determine the
#' quality of an item. It is defined as
#' \deqn{S = \frac{r_{it}}{2 \times{} \sqrt{Difficulty \times{} (1-Difficulty)}}}
#'
#' @note \code{klausur} allows some tweaks that are probably not as useful as they seem. For instance, having items with more than one correct answer doesn't
#' necessarily yield more diagnostic information, allowing for those being answered partially adds to that, and binding marks blindly to a normal distribution can
#' give quite unfair test results! In addition, please do \strong{always check a sample of the results} to make sure no errors accurred.
#' 
#' @param data An object of class \code{\link[klausuR]{klausuR.answ-class}}.
#' @param marks A vector assigning marks to points achieved (see details). Alternatively, set it to \code{"suggest"} to let
#'    \code{\link[klausuR:klausur.gen.marks]{klausur.gen.marks}} calculate suggestions under the assumption of normal distribution.
#'    If \code{NULL}, this value must be set in the \code{data} object.
#' @param mark.labels If \code{marks="suggest"}, use these as the marks you want to give.
#' @param items Indices of a subset of variables in \code{data} to be taken as items.
#' @param wght A vector with weights for each item (named also according to \code{Item###}). If \code{NULL}, the value from the \code{data} object
#'    will be used.
#' @param score Specify the scoring policy, must be one of \code{"solved"} (default), \code{"partial"}, \code{"liberal"}, \code{"pick-n"},
#'    \code{"NR"}, \code{"ET"}, \code{"NRET"}, or \code{"NRET+"}.
#' @param matn A matriculation number of a subject, to receive detailed results for that subject.
#' @param na.rm Logical, whether cases with NAs should be ignored in \code{data}. Defaults to TRUE.
#' @param cronbach Logical. If TRUE, Cronbach's alpha will be calculated.
#' @param item.analysis Logical. If TRUE, some usual item statistics like difficulty, discriminatory power and distractor analysis will be calculated.
#' @param sort.by A character string naming the variable to sort the results by. Set to \code{c()} to skip any re-ordering.
#'  If \code{cronbach} is TRUE, too, it will include the alpha values if each item was deleted.
#' @param maxp Optional numeric value, if set will be forced as the maximum number of points achievable. This should actually not be needed,
#'    if your test has no strange errors. But if for example it later turns out you need to adjust one item because it has two instead of
#'    one correct answers, this option can become handy in combination with "partial" scoring and item weights.
#' @return An object of class \code{\link[klausuR]{klausuR-class}} with the following slots.
#'  \item{results}{A data.frame with global results}
#'  \item{answ}{A data.frame with all given answers}
#'  \item{corr}{A vector with the correct answers}
#'  \item{wght}{A vector with the weights of items}
#'  \item{points}{A data.frame with resulting points given for the answers}
#'  \item{marks}{A vector with assignments of marks to achieved score}
#'  \item{marks.sum}{A more convenient matrix with summary information on the defined marks}
#'  \item{trfls}{A data.frame of TRUE/FALSE values, whether a subject was able to solve an item or not}
#'  \item{anon}{A data.frame for anonymous feedback}
#'  \item{mean}{A table with mean, median and quartiles of the test results}
#'  \item{sd}{Standard deviation of the test results}
#'  \item{cronbach}{Internal consistency, a list of three elements "alpha", "ci" (confidence interval 95\%) and "deleted" (alpha if item was removed)}
#'  \item{item.analysis}{A data.frame with information on difficulty, discriminatory power, discriminant factor and Lienert's selection index of all items.}
#'  \item{distractor.analysis}{A list with information on the selected answer alternatives for each individual item (only calculated if \code{item.analysis=TRUE}).
#'    Also lists the discriminatory power of each alternative, being the point-biserial (a.k.a Pearson) correlation of it with the global outcome.}
#'  \item{misc}{Anything that was stored in the \code{misc} slot of the input data.}
#'  Not all slots are shown by default (refer to \code{\link[klausuR:show-methods]{show}} and \code{\link[klausuR:plot]{plot}}).
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @seealso \code{\link[klausuR:klausur.data]{klausur.data}}, \code{\link[klausuR:klausur.report]{klausur.report}},
#'    \code{\link[klausuR:compare]{compare}},
#'    \code{\link[klausuR:klausur.gen]{klausur.gen}}, \code{\link[klausuR:klausur.gen.marks]{klausur.gen.marks}},
#'    \code{\link[klausuR:klausur.gen.corr]{klausur.gen.corr}},
#'    \code{\link[klausuR:plot]{plot}}, \code{\link[psychometric]{psychometric}}
#' @keywords misc
#' @import psychometric
#' @import polycor
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
#' 
#' # to try pick-n scoring, we must also define all distractors
#' falsch <- c(Item01=1245, Item02=1345, Item03=1345, Item04=1345, Item05=1235,
#'  Item06=1245, Item07=1235, Item08=2345, Item09=1345, Item10=1345, Item11=1235,
#'  Item12=1235, Item13=1345, Item14=1245, Item15=1345, Item16=1245, Item17=1235,
#'  Item18=1235, Item19=1245, Item20=1234, Item21=1245, Item22=1245, Item23=2345,
#'  Item24=1245, Item25=2345, Item26=1245, Item27=1234, Item28=1245, Item29=1235,
#'  Item30=1235, Item31=245, Item32=15)
#'
#' data.obj <- klausur.data(answ=antworten, corr=richtig, wrong=falsch, marks=notenschluessel)
#' klsr.obj <- klausur(data.obj, score="pick-n")
#'
#' ############################
#'  # example for an NRET test
#' ############################
#' # load sampla data in SPSS format
#' data(spss.data)
#' # define correct answers
#' spss.corr <- c(
#'    item01=2, item02=3, item03=3, item04=3, item05=2,
#'    item06=2, item07=3, item08=1, item09=1, item10=2)
#'
#' # convert into klausuR type coding
#' klausuR.data <- nret.translator(spss.data, spss="in")
#' klausuR.corr <- nret.translator(spss.corr, spss="in", corr=TRUE,
#'   num.alt=3, spss.prefix=c(corr="item"))
#' # now create the data object; "Nickname" must be renamed
#' data.obj <- klausur.data(answ=klausuR.data, corr=klausuR.corr,
#'   rename=c(Pseudonym="Nickname"))
#' 
#'  # finally, the test can be evaluated, using the scoring functions available
#' NRET.results <- klausur(data.obj, marks="suggest", mark.labels=11, score="NRET")
#' NRETplus.results <- klausur(data.obj, marks="suggest", mark.labels=11, score="NRET+")
#' NR.results <- klausur(data.obj, marks="suggest", mark.labels=11, score="NR")
#' ET.results <- klausur(data.obj, marks="suggest", mark.labels=11, score="ET")

klausur <- function(data, marks=NULL, mark.labels=NULL, items=NULL, wght=NULL, score="solved",
    matn=NULL, na.rm=TRUE, cronbach=TRUE, item.analysis=TRUE, sort.by="Name", maxp=NULL){
    ## whenever these options/arguments should change, update klausur.mufo() accordingly!

    if(!inherits(data, "klausuR.answ")){
      stop(simpleError("'data' must be of class 'klausuR.answ'!"))
    } else {}

    corr <- slot(data, "corr")[["corr"]]
    wrong <- slot(data, "corr")[["wrong"]]
    if(!is.null(slot(data, "corr")[["corr.key"]])){
      warning("This test seems to have multiple test forms. Perhaps try klausur.mufo() instead?", call.=FALSE)
    } else {}

    if(is.null(wght)){
      wght <- slot(data, "score")[["wght"]]
    } else {}

    if(is.null(marks)){
      if(is.null(slot(data, "score")[["marks"]])){
        stop(simpleError("You must give some value for 'marks', either in 'data' or with the klausur() call!"))
      } else {
        marks <- slot(data, "score")[["marks"]]
      }
    } else {}

    ## TODO: clean up checks, so cbind is not needed!
    ## firstly, check input data and quit if necessary
    # data.check.klausur() is an internal function, defined in klausuR-internal.R
    sane.data <- data.check.klausur(answ=cbind(slot(data, "id"), slot(data, "items")), corr=corr, items=items, na.rm=na.rm)
    stopifnot(scoring.check.klausur(corr=corr, marks=marks, wght=wght, score=score, maxp=maxp))
    answ <- sane.data$answ
    items <- sane.data$items

    # probably sort cases
    if(length(sort.by) > 0){
      if(!sort.by %in% names(answ)){
        stop(simpleError(paste("Can't sort by '",sort.by,"', there's no such variable!", sep="")))
      } else {}
      new.order <- order(answ[[sort.by]])
      answ <- answ[new.order,]
    } else {}

    ### results section
    # set default min.score to zero
    min.score <- 0
    # probably weight items, to calculate the maximum score
    # NRET et al. can't be weighted yet
    if(is.null(wght) | score %in% c("NR", "ET", "NRET", "NRET+")){
      nret.test.chars <- nret.minmax(corr=corr, score=score)
      if(is.null(maxp)){
        if(is.null(data@score$maxp)){
          maxp <- nret.test.chars["maxp"]
        } else {
          maxp <- data@score$maxp
          warning(paste("Forced maximum number of points to ", maxp, " (as set in data object)!", sep=""), call.=FALSE)
        }
      } else {
        warning(paste("Manually forced maximum number of points to ", maxp, "!", sep=""), call.=FALSE)
      }
      min.score <- nret.test.chars["minp"]
      baseline <- nret.test.chars["baseline"]
      num.alt <- nret.test.chars["num.alt"]
      # for the results, create a vector of 1's by number of items
      wght.results <- rep(1, length(items))
    } else {
      if(is.null(maxp)){
        if(is.null(data@score$maxp)){
          maxp <- sum(wght)
        } else {
          maxp <- data@score$maxp
          warning(paste("Forced maximum number of points to ", maxp, " (as set in data object)!", sep=""), call.=FALSE)
        }
      } else {
        warning(paste("Manually forced maximum number of points to ", maxp, "!", sep=""), call.=FALSE)
      }
      wght.results <- wght
    }

    # create the TRUE/FALSE-matrix for solved items
    # this will become a point matrix instead if partial results are ok
    if(length(score) == 1 && score %in% c("partial", "pick-n", "liberal", "NR", "ET", "NRET", "NRET+")){
      # weights will be considered here as well
      # if partial answers should be considered
      if(identical(score, "partial")){
        wahr.falsch <- data.frame(sapply(items, function(x){
            return(partial(item.answ=answ[x], corr=corr, wght=1, strict=TRUE, mode="percent"))
          }), stringsAsFactors=FALSE)
        ergebnisse  <- data.frame(sapply(items, function(x){
            return(partial(item.answ=answ[x], corr=corr, wght=wght[which(items == x)], strict=TRUE, mode="percent"))
          }), stringsAsFactors=FALSE)
      } else if(identical(score, "pick-n")){
        wahr.falsch <- data.frame(sapply(items, function(x){
            return(pickN(item.answ=answ[x], corr=corr, wrong=wrong, wght=1, mode="percent"))
          }), stringsAsFactors=FALSE)
        ergebnisse  <- data.frame(sapply(items, function(x){
            return(pickN(item.answ=answ[x], corr=corr, wrong=wrong, wght=wght[which(items == x)], mode="percent"))
          }), stringsAsFactors=FALSE)
      } else if(identical(score, "liberal")){
        wahr.falsch <- data.frame(sapply(items, function(x){
            return(partial(item.answ=answ[x], corr=corr, wght=1, strict=FALSE, mode="percent"))
          }), stringsAsFactors=FALSE)
        ergebnisse  <- data.frame(sapply(items, function(x){
            return(partial(item.answ=answ[x], corr=corr, wght=wght[which(items == x)], strict=FALSE, mode="percent"))
          }), stringsAsFactors=FALSE)
      } else if(score %in% c("NR", "ET", "NRET", "NRET+")){
        wahr.falsch <- data.frame(sapply(items, function(x){
            return(nret.score(answ[[x]], corr=corr[names(answ[x])], score=score, is.true="+", is.false="-", missing="0", err="*",
            num.alt=num.alt, true.false=TRUE))
          }), stringsAsFactors=FALSE)
        ergebnisse  <- data.frame(sapply(items, function(x){
            return(nret.score(answ[[x]], corr=corr[names(answ[x])], score=score, is.true="+", is.false="-", missing="0", err="*",
            num.alt=num.alt, true.false=FALSE))
          }), stringsAsFactors=FALSE)
      } else {}
      colnames(wahr.falsch) <- colnames(answ[items])
      colnames(ergebnisse)  <- colnames(answ[items])
    } else {
      wahr.falsch <- data.frame(t(t(answ[,items]) == corr))
      # in case weights were defined, e.g. for items with multiple correct answers, take them into account
      if(!is.null(wght)) {
        ergebnisse <- data.frame(t(t(wahr.falsch) * wght))
      } else {
        # we'll add 0 to forcibly convert logical values to numerics
        ergebnisse <- wahr.falsch + 0
      }
    }
    # make a copy of the matrix with resulting points, especially interesting if items were weighted
    ergebnisse.daten <- cbind(MatrNo=answ$MatrNo, ergebnisse)

    ## calculate points
    punkte <- rowSums(ergebnisse)

    ## descriptive statistics
    mittel.quart <- summary(punkte)
    stdabw <- sd(punkte, na.rm=TRUE)

    # should marks be suggested or are they given?
    if(length(marks) == 1 && identical(marks, "suggest")){
      warning("Marks are a suggestion, not your definition!", call.=FALSE)
      ##
      if(score %in% c("ET", "NRET", "NRET+")){
        # a hack to force the right max. points here
        marks <- klausur.gen.marks(mark.labels=mark.labels, answ=maxp, wght=NULL, suggest=list(mean=mean(punkte), sd=stdabw), minp=baseline)
        } else if(identical(score, "NR")){
          marks <- klausur.gen.marks(mark.labels=mark.labels, answ=maxp, wght=NULL, suggest=list(mean=mean(punkte), sd=stdabw), minp=min.score)
      } else {
        marks <- klausur.gen.marks(mark.labels=mark.labels, answ=answ, wght=wght, suggest=list(mean=mean(punkte), sd=stdabw), minp=min.score)
      }
    }

    # check wheter maximum score matches the assigned marks
    if(length(marks) != maxp){
      warning(paste("Achievable score and marks do not match!\n  Maximum score of test:",maxp,"\n  Maximum score assigned to marks:",length(marks)), call.=FALSE)
    }

    # assign the marks. zero points will result in lowest mark, NA will stay NA
    note <- sapply(punkte, function(x){
        if(is.na(x))
          return(NA) else{}
        if(x > 0)
          note <- marks[x]
        else
          note <- marks[1]
        return(note)
      })

    # create data object with name, mat-nr and global results
    # calls the internal function global.results()
    if(score %in% c("ET", "NRET", "NRET+")){
      ergebnis.daten <- global.results(answ=answ, points=punkte, maxp=maxp, mark=note, minp=baseline)
    } else {
      ergebnis.daten <- global.results(answ=answ, points=punkte, maxp=maxp, mark=note, minp=min.score)
    }
    # if pseudonyms were given, use them for anonymous feedback
    ergebnis.anonym <- anon.results(ergebnis.daten)

    ## psychometic quality of the items
    if(score %in% c("ET", "NRET", "NRET+")){
      # if ET/NRET type items are to be analyzed, try with points instead of true/false matrix
      # this has implication for interpretation
      item.values.to.anl <- ergebnisse
      if(isTRUE(item.analysis)){
        warning("Be careful with interpreting item analysis on ET/NRET data.\n  It's not dichotomously scaled, analysis was done with the achieved scores.", call.=FALSE)
      } else {}
    } else {
      item.values.to.anl <- wahr.falsch
    }
      if(isTRUE(cronbach)){
        # calling an internal function which is
        # using alpha() from package "psychometric"
        cron.alpha.list <- calc.cronbach.alpha(na.omit(item.values.to.anl))
      } else {
        cron.alpha.list <- list(alpha=NULL, ci=NULL, deleted=NULL)
      }
      if(isTRUE(item.analysis)){
        # calling another internal function which is also
        # using alpha() from package "psychometric"
        item.analyse <- calc.item.analysis(item.values.to.anl, cron.alpha.list)
        distractor.analysis <- distrct.analysis(
            answ=slot(data, "items"),
            corr=corr,
            points=ergebnisse.daten,
            results=data.frame(MatrNo=answ$MatrNo, Points=ergebnis.daten$Points),
            partWhole=isTRUE(score %in% c("solved", "NR"))
          )
      } else {
        item.analyse <- data.frame(NULL)
        distractor.analysis <- list(NULL)
      }

    ## compose the resulting object
    # here we make a copy of the TRUE/FALSE matrix, to be able to check each result individually
    wahrfalsch.daten <- cbind(MatrNo=answ$MatrNo, wahr.falsch)
    # the same way make a copy with given answers as well
    antwort.daten <- cbind(MatrNo=answ$MatrNo, answ[,items])
    # use the internal marks.summary() function to create convenient information on the mark definitions
    marks.info <- marks.summary(marks, minp=min.score)
    # here we finally take all the partial results an knit them together into one object
    alle.ergebnisse <- new("klausuR",
          results=ergebnis.daten,
          answ=antwort.daten,
          corr=corr,
          wght=wght.results,
          points=ergebnisse.daten,
          marks=marks,
          marks.sum=marks.info,
          trfls=wahrfalsch.daten,
          anon=ergebnis.anonym,
          mean=mittel.quart,
          sd=stdabw,
          cronbach=cron.alpha.list,
          item.analysis=item.analyse,
          distractor.analysis=distractor.analysis,
          misc=data@misc)

    ## output options
    # return all data or just for one matriculation number?
    if(!is.null(matn)){
      if(!sum(answ$MatrNo == matn) == 1){
      stop(simpleError("The given matriculation number is not defined!"))
      }
      pers.ergebnisse <- ergebnis.daten[ergebnis.daten$MatrNo == matn,]
      items.solved <- wahrfalsch.daten[wahrfalsch.daten$MatrNo == matn, -1]
      rownames(items.solved) <- "Solved"
      ausgabe <- list(Results=pers.ergebnisse, SolvedItems=t(items.solved))
    }
    else
        # if matn is not set, return the default result object
        ausgabe <- alle.ergebnisse

    ## ausgabe
    return(ausgabe)
    }
