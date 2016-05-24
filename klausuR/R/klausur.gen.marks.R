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


#' Generate mark assignments
#'
#' Create a vector of marks to be used by \code{\link[klausuR:klausur]{klausur}}.
#'
#' If \code{mark.labels} is set to one of the arguments 6, 11, 16, "DIHK", "UK", "USA" or "A", often used schemes for marks will be used as
#' a preset. In case of \code{mark.labels=6}, marks will go from 1 (best) to 6 (failed), as widely used in German schools.
#' In case of \code{mark.labels=11}, marks will range from 1.0 (best) to 5.0 (failed), with the marks 1 through 3 being
#' split into three decimal steps .0, .3 and .7, as is often used in academic institutions. In case of \code{mark.labels=16},
#' marks will be a range of points from 15 (best) to 0 (failed), as often used in German gymnasiums. If \code{mark.labels="A"},
#' marks A to F are given.
#'
#' For the other cases some more probably useful assumptions are being made, which percentage of achieved points leads to which mark.
#' If \code{mark.labels="DIHK"}, marks will be 1 through 6, and  calculated according to usual standards of the Deutsche Industrie- und Handelskammer
#' (1 > 92\%, 2 > 81\%, 3 > 67\%, 4 > 50\%, 5 > 30\%, 6 below that).
#' If \code{mark.labels="UK"}, marks are A > 90\%, B > 65\%, C > 35\%, D > 10\% and E below that, and for \code{mark.labels="USA"} it's
#' A > 90\%, B > 80\%, C > 70\%, D > 60\% and F below that. Please note that the percentages indicate individual test results
#' and not "the best X percent of the sample".
#' If you'd rather use your own system, either declare it as a vector, or leave as NULL, and you'll be
#' asked (be sure to \strong{begin with the worst} mark!).
#'
#' The parameter \code{answ} is quite versatile as well. You can just feed it your observation data, if it complies
#' with the naming scheme for items (\code{Item###}, see also \code{\link[klausuR:klausur.gen]{klausur.gen}}), and
#' \code{klausur.gen.marks} will calculate the maximum score automatically. Or you assign the maximum directly as an integer
#' value.
#'
#' Another feature can be toggled with the parameter \code{suggest}. If you feed it with the mean and standard deviation values of your test's
#' results, marks are automatically assigned to the achieved score under the assumption of normal distribution. Please understand that
#' the naming "suggest" is not an accident! This is only a suggestion, please review it, tweak it, revise it, until it fits your needs.
#' However, this feature can directly be called by \code{\link[klausuR:klausur]{klausur}}.
#'
#' @param mark.labels Either a vector with labels (names) for all marks that should be assigned to certain test scores,
#'        or one of the arguments 6, 11, 16, "DIHK", "UK", "USA" or "A" (see Details). If NULL, you will be asked to type in labels.
#' @param answ Either an object with item names in klausuR scheme (see \code{\link[klausuR:klausur]{klausur}}),
#'        e.g. your observation data, or an integer representing the maximum score of the test. If NULL, you will
#'        be asked for the maximum score.
#' @param wght A vector with weights for each item. If NULL, the number of items is used as maximum score, that is,
#'        each item gives a point.
#' @param suggest A list with the elements \code{mean} and \code{sd}. If both are not NULL, this function will suggest marks for achieved points
#'    assuming normal distribution. That is, "mean" and "sd" should be set to the corresponding values of the test's results.
#' @param minp An integer value, in case there is a minimum score no-one can fall below (which can happen, e.g., with ET/NRET scoring and
#'    different numbers of answer alternatives). Should be left as is in most cases.
#' @return A character vector.
#' @keywords utilities
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @examples
#' \dontrun{
#' notenschluessel <- klausur.gen.marks(mark.labels=11,answ=antworten)
#' }
#' @seealso \code{\link[klausuR:klausur]{klausur}}
#' @export

klausur.gen.marks <- function(mark.labels=NULL, answ=NULL, wght=NULL, suggest=list(mean=NULL, sd=NULL), minp=0){

  ## function max.score
  # used to determine the maximum score
  max.score <- function(answ, wght){
    # if weights are given, use them to calculate maximum score:
    if(!is.null(wght)){
      maxp <- sum(wght)
    } else {
      # otherwise use the number of items
      if(!is.null(answ)){
        if(is.double(answ) && length(answ) == 1){
          maxp <- answ
        } else {
          maxp <- as.numeric(length(grep("Item([[:digit:]]{1,3})",names(answ))))
        }
      } else {
        maxp <- as.numeric(readline(paste("What is the maximum score achievable?")))
      }
    }
    # stop if still not useful
    if(length(maxp) != 1){
      stop(simpleError(paste("Illegal value for maximum score: \"",maxp,"\"", sep="")))
    } else {
      cat("Maximum score:\t", maxp,"\n")
      return(maxp)
    }
  } ## end function max.score

  ## function label.marks
  # setting label marks, what a surprise...
  label.marks <- function(mark.labels){
    if(is.null(mark.labels)){
      mark.labels <- c()
      markname <- "placeholder"
      mark.labels[1] <- as.character(readline(paste("Please start with the worst mark you'll assign:")))
      index <- 2
      while(!identical(markname, "")){
        markname <- as.character(readline(paste("Please name the next better mark (leave empty to finish!):")))
        if(!identical(markname, "")){
          mark.labels[index] <- markname
          index <- index + 1
        } else{}
      }
    } else {
      # check if a preset of marks can be used
      if(length(mark.labels) == 1 && identical(as.numeric(mark.labels), 16)){
        mark.labels <- c("00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15")
      } else if(length(mark.labels) == 1 && identical(as.numeric(mark.labels), 11)){
        mark.labels <- c("5.0","4.0","3.7","3.3","3.0","2.7","2.3","2.0","1.7","1.3","1.0")
      } else if(length(mark.labels) == 1 && (identical(as.numeric(mark.labels), 6) || identical(mark.labels, "DIHK"))){
        mark.labels <- c(6:1)
      } else if(length(mark.labels) == 1 && identical(mark.labels, "USA")){
        mark.labels <- c("F","D","C","B","A")
      } else if(length(mark.labels) == 1 && identical(mark.labels, "UK")){
        mark.labels <- c("E","D","C","B","A")
      } else if(length(mark.labels) == 1 && identical(mark.labels, "A")) {
        mark.labels <- c("F","E","D","C","B","A")
      } else{}
    }
    # stop if still not useful
    if(is.null(mark.labels) || length(mark.labels) == 1){
      stop(simpleError(paste("Illegal value for mark labels: \"",mark.labels,"\"", sep="")))
    } else {
      cat("Mark labels:\t", mark.labels,"\n")
      return(mark.labels)
    }
  } ## end function label.marks

  ## function read.points
  # this function is called to ask for the assigned points
  read.points <- function(markname, maxp){
    points <- ""
    while(identical(points, "")){
      points <- readline(paste("Please input the minimum score for mark ",as.character(markname),":", sep=""))
      if(!identical(points, "")){
        if(points <= maxp){
          if(points > points.assigned) {
            marks[points:maxp] <<- as.character(markname)
            cat("OK: Assigned mark \"",markname,"\" to a minimum score of ",points," points.\n", sep="")
            points.assigned <<- points
          } else {
            cat("Illegal: Point value (",points,") already assigned to mark \"",marks[as.numeric(points)],"\"!\n", sep="")
            points <- ""
          }
        } else {
          cat("Illegal: Point value (",points,") exceeding maximum score!\n", sep="")
          points <- ""
        }
      } else{}
    }
  } ## end function read.points

    ## function suggestion
    # this function takes information un mean and standard deviation of test results,
    # and suggests marks according to normal distribution
    suggestion <- function(mark.labels, maxp, mean, sd, minp){
      # so, which quantiles do we need? obviously, one less that we have marks
      # we take one more than the number of marks, because the first and last will be 0 and 1,
      # so we drop the last one again and keep the rest
      quants <- seq(from=0, to=1, by=(1/length(mark.labels))) [c(-(length(mark.labels)+1))]

      # now have a look at the normal distribution,
      # get the quantiles for mean and sd of our test results
      lapply(mark.labels[-1], function(x){
        # using max() to make sure min.points can't fall below 1!
        min.points <- max(max(1, minp), ceiling(qnorm(quants[which(mark.labels == x)], mean=mean, sd=sd)))
        # in case results are odd, let's at least not have unrealistic points values:
        if(min.points > maxp){
          min.points <- maxp
        } else {}
          marks[min.points:maxp] <<- x
        }
      )
      return(marks)
    } ## end function suggestion

  ## function scheme.calc
  scheme.calc <- function(scheme, maxp, mark.labels.set, marks){
    if(identical(scheme, "DIHK")){
      quants <- c(0, .30, .50, .67, .81, .92)
    } else if(identical(scheme, "USA")){
      quants <- c(0, .60, .70, .80, .90)
    } else if(identical(scheme, "UK")){
      quants <- c(0, .10, .35, .65, .90)
    } else {
      return(NA)
    }

    # calculate threshold values, beginning with 1 point, ending with maximum
    mark.thresholds <- ceiling(maxp * quants)
    mark.thresholds[1] <- 1
    mark.thresholds <- append(mark.thresholds, maxp)
    # fill the mark vector
    marks <- lapply(c(1:length(quants)), function(x){
      lower.lim <- mark.thresholds[x]
      upper.lim <- mark.thresholds[x+1]
      marks[lower.lim:upper.lim] <<- mark.labels.set[x]
      return(marks)
    })

    return(unlist(marks[1]))
  }  ## end function scheme.calc


  ## now let's do it!!!
  # call max.score() to calculate maximum score
  maxp <- max.score(answ, wght)
  # rename the marks
  mark.labels.set <- label.marks(mark.labels)

  # create empty object for the assignments
  marks <- c()
  # and fill it with the worst mark for a start
  marks[1:maxp] <- as.character(mark.labels.set[1])

  # now read it all in!
  if(identical(mark.labels, "DIHK") || identical(mark.labels, "USA") || identical(mark.labels, "UK")){
    marks <- scheme.calc(mark.labels, maxp, mark.labels.set, marks)
  } else if(!is.numeric(suggest[["mean"]]) || !is.numeric(suggest[["sd"]])) {
    # set counter to zero
    points.assigned <- 0
    # then call the function that uses it
    mapply(read.points, mark.labels.set[-1], MoreArgs=list(maxp=maxp))
  } else {
    marks <- suggestion(mark.labels.set, maxp, mean=suggest[["mean"]], sd=suggest[["sd"]], minp=minp)
  }

  return(marks)
}
