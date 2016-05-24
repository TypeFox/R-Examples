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


#' Plot methods for S4 objects of class klausuR and klausuR.mult
#'
#' These plot methods are beeing called by \code{\link[klausuR:klausur.report]{klausur.report}}.
#' If \code{x} is of class \code{klausuR.mult}, only the global results will be plotted.
#' Should you rather like plots on each test form, call \code{plot} with the single slots from that object accordingly.
#'
#' @aliases plot,-methods plot,klausuR-method plot,klausuR,missing-method plot,klausuR.mult-method plot,klausuR.mult,missing-method
#' @param x An S4 object of class \code{klausuR} or \code{klausuR.mult}
#' @param y From the generic \code{plot} function, ignored for klausuR class objects.
#' @param marks Logical, whether the histogram should show the distribution of points (default) or marks
#' @param sd.lines Logical, whether standard deviation lines should be plotted
#' @param plot.normal Logical, whether normal distribution should be plotted (according to mean and Sd of the results)
#' @param na.rm Logical, whether NA values should be ignored. Defaults to TRUE, because plotting would fail otherwise
#' @param ... Any other argument suitable for plot()
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @seealso \code{\link[klausuR:klausur]{klausur}}, \code{\link[klausuR:klausur.mufo]{klausur.mufo}}, \code{\link[klausuR:klausur.report]{klausur.report}}
#' @keywords methods plot
#' @import methods
#' @export
#' @docType methods
#' @rdname plot-methods
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
#' data.obj <- klausur.data(answ=antworten, corr=richtig, marks=notenschluessel)
#' klsr.obj <- klausur(data.obj)
#' plot(klsr.obj, marks=TRUE)
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' @docType methods
#' @aliases plot,klausuR,missing,ANY-method
#' @rdname plot-methods
setMethod("plot", signature(x="klausuR", y="missing"), function(x, marks=FALSE, sd.lines=FALSE, plot.normal=TRUE, na.rm=TRUE, ...){

  klsr <- x
  if(isTRUE(na.rm)){
    erg.points <- na.omit(klsr@results$Points)
    erg.marks <- na.omit(klsr@results$Mark)
  } else {
    erg.points <- klsr@results$Points
    erg.marks <- klsr@results$Mark
  }

  if(isTRUE(marks)){
    erg.max.freq <- max(summary(as.factor(erg.marks)))
    erg.marks.num <- erg.marks
    # let's check if marks are a character factor
    if(is.character(erg.marks)){
      erg.mark.levels <- levels(as.factor(erg.marks))
      # now we'll replace them with numbers, that is we assume an ordinal scale
      lapply(erg.mark.levels, function(x){erg.marks.num[erg.marks.num == x] <<- which(erg.mark.levels == x)})
      erg.marks.num <- as.numeric(erg.marks.num)
    } else {}
    norm.mean <- mean(erg.marks.num)
    norm.sd <- sd(erg.marks.num)
    norm.min <- min(erg.marks.num)
    norm.max <- max(erg.marks.num)

    plot(as.ordered(erg.marks),yaxp=c(0,erg.max.freq,erg.max.freq),...)
    if(plot.normal){
      # plot normal distribution
      par(new=TRUE)
      plot(function(x) dnorm(x, mean=norm.mean, sd=norm.sd), from=norm.min, to=norm.max, axes=FALSE, xlab="", ylab="", lwd=5, lty=3)
    } else {}
  } else {
    erg.min <- min(erg.points)
    erg.max <- max(erg.points)
    erg.mittel <- mean(erg.points)
    erg.sd <- sd(erg.points)
    hist.breaks <- axis.breaks(erg.min, erg.max)
    erg.max.freq <- max(summary(as.factor(round(erg.points))))
    
    hist(erg.points, breaks=hist.breaks[["breaks"]], col="grey", xaxt="n", yaxp=c(0,erg.max.freq,erg.max.freq), ...)
    axis(1, at=c((hist.breaks[["min"]]+0.5):(hist.breaks[["max"]]-0.5)), labels=c((hist.breaks[["min"]]+1):hist.breaks[["max"]]))
    if(plot.normal){
      # plot normal distribution
      par(new=TRUE)
      plot(function(x) dnorm(x, mean=erg.mittel, sd=erg.sd), from=erg.min, to=erg.max, axes=FALSE, xlab="", ylab="", lwd=5, lty=3)
    } else {}
    # if desired, plot standard deviation
    if(sd.lines){
      abline(v=c(erg.mittel-2*erg.sd, erg.mittel-erg.sd, erg.mittel, erg.mittel+erg.sd, erg.mittel+2*erg.sd), lwd=1, lty=1)
    } else {}
  }
})

#' @docType methods
#' @aliases plot,klausuR.mult,missing,ANY-method
#' @rdname plot-methods
setMethod("plot", signature(x="klausuR.mult", y="missing"), function(x, marks=FALSE, sd.lines=FALSE, plot.normal=TRUE, ...){
  klausur.global.object <- x@results.glob
  plot(x=klausur.global.object, marks=marks, sd.lines=sd.lines, plot.normal=plot.normal, ...)
})
