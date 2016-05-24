##' A \code{print} method for \code{xBalance} objects.
##'
##' @title Printing xBalance Objects
##' @param x An object of class "xbal" which is the result of a call
##'   to \code{xBalance}.
##' @param which.strata The stratification candidates to include in
##'   the printout. Default is all.
##' @param which.stats The test statistics to include. Default is all
##'   those requested from the call to \code{xBalance}.
##' @param which.vars The variables for which test information should
##'   be displayed. Default is all.
##' @param print.overall Should the omnibus test be reported? Default
##'   is \code{TRUE}.
##' @param digits To how many digits should the results be displayed?
##'   Default is \code{max(2,getOptions("digits")-4)}.
##' @param printme Print the table to the console? Default is
##'   \code{TRUE}.
##' @param show.signif.stars Use stars to indicate z-statistics larger
##'   than conventional thresholds. Default is \code{TRUE}.
##' @param show.pvals Instead of stars, use p-values to summarize the
##'   information in the z-statistics. Default is \code{FALSE}.
##' @param horizontal Display the results for different candidate
##'   stratifications side-by-side (Default, \code{TRUE}), or as a
##'   list for each stratification (\code{FALSE}).
##' @param ... Other arguements. Not currently used.
##' @return \describe{
##' \item{vartable}{The formatted table of variable-by-variable
##'   statistics for each stratification.}
##' \item{overalltable}{If the overall Chi-squared statistic is
##'   requested, a formatted version of that table is returned.}
##' }
##' @seealso xBalance
##' @export
##' @aliases print
##' @keywords print
##' @examples
##' data(nuclearplants)
##'
##' xb0<-xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'               data=nuclearplants)
##'
##' print(xb0)
##'
##' xb1<-xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          strata = data.frame(unstrat = factor(character(32)),
##'                              pt = factor(nuclearplants$pt)),
##'          data = nuclearplants,
##'          report = c("all"))
##'
##' str(xb1)
##'
##' print(xb1)
##'
##' print(xb1, show.pvals = TRUE)
##'
##' print(xb1, horizontal = FALSE)
##'
##' ## The following doesn't work yet.
##' \dontrun{print(xb1, which.vars=c("date","t1"),
##'          which.stats=c("adj.means","z.scores","p.values"))}
##'
##' ## The following example prints the adjusted means
##' ## labeled as "treatmentvar=0" and "treatmentvar=1" using the
##' ## formula provided to xBalance().
##'
##' print(xb1,
##'       which.vars = c("date", "t1"),
##'       which.stats = c("pr=0", "pr=1", "z", "p"))
##'
##' ## Now, not asking for the omnibus test
##' xb2 <- xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          strata = data.frame(unstrat = factor(character(32)),
##'                              pt = factor(nuclearplants$pt)),
##'          data = nuclearplants,
##'          report = c("all"))
##'
##' print(xb2, which.strata = "pt")
print.xbal <- function (x, which.strata=dimnames(x$results)[["strata"]],
                        which.stats=dimnames(x$results)[["stat"]],
                        which.vars=dimnames(x$results)[["vars"]],
                        print.overall=TRUE,
                        digits = NULL, printme=TRUE,
                        show.signif.stars=getOption("show.signif.stars"),
                        show.pvals=!show.signif.stars,
                        horizontal=TRUE,...) {
  ##Notes: right now we've decided that you can't print both signif stars and p-values. make a choice.

  # withOptions will allow us to safely reset the digits
  # even if an error is thrown the option should be the same after this function
  DIGITS = ifelse(is.null(digits), max(2, getOption("digits")-4), digits)
  withOptions(list(digits = DIGITS), function() {
    ##makeSigStarsStdNormal <- function(zs) {
    ##  if (length(zs)){##c('','.  ','*  ','** ','***'
    ##    factor(c('','.','*','**','***') [1+
    ##             apply(abs(zs)>=matrix(qnorm(c(.95, .975, .995,.9995)),
    ##                        length(zs),4, byrow=TRUE),1, sum)]
    ##           )} else {character(0)}
    ##                                     }

    theresults <- x$results
    if("overall" %in% names(x)) { ##Extract the omnibus chisquared test from the xbal object...
      theoverall <- x$overall
    } else {
      theoverall<-NULL ##..or set it to NULL if it does not exist.
    }

    if (length(theresults) == 0 & is.null(theoverall)) {
      stop("There is a problem. Probably all of the variables (",
           all.vars(formula(x)),
           ") are constants within strata. Or else there is some other problem, try debug(RItools:::xBalance) to see what might be going on.")
    }

    if (length(theresults)==0 & !is.null(theoverall)){##The user has requested only the omnibus test and not the tests for the individual variables
      theresults<-NULL
      thevartab<-NULL
    }

    signifier <- function(data) {
      symnum(data,
             corr = FALSE,
             na = FALSE,
             abbr.colnames=FALSE, ##from print.summary.lm
             cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
             symbols = c("***", "** ", "*  ", ".  ", "   "))
    }

    ftabler <- function(data) {   ##Summarize the variable-by-variable output array as a flat contingency table
      ftable(data, col.vars=c("strata","stat"),row.vars=c("vars"))
    }

    if (show.signif.stars & !show.pvals & !is.null(theresults)) {

      Signif <- signifier(theresults[,"p",which.strata,drop=FALSE])

      ##Nicer alignment, but not as nice labels
      ##junk<-do.call(cbind,lapply(which.strata,function(x){cbind(as.data.frame(theresults[,,x])," "=format(Signif[,,x]))}))
      newresults <- array(dim=dim(theresults)+c(0,1,0),
                          dimnames=list(vars=dimnames(x$results)[["vars"]],
                              stat=c(dimnames(x$results)[["stat"]],"sig."),
                              strata=dimnames(x$results)[["strata"]]))

      newresults[,-grep("sig.", dimnames(newresults)[[2]]),] <- format(theresults,DIGITS)
      newresults[dimnames(Signif)[["vars"]], "sig.",dimnames(Signif)[["strata"]]]<-format(Signif)


      if (horizontal){
        theftab <- ftabler(
            newresults[which.vars, c(which.stats[!(which.stats=="p")],"sig."), which.strata, drop=FALSE])

        attr(theftab,"col.vars")$stat[attr(theftab,"col.vars")$stat=="sig."] <- ""
        if ("z" %in% which.stats) {
          attr(theftab,"col.vars")$stat[attr(theftab,"col.vars")$stat=="z"]<-"   z   "
        }

        thevartab<-theftab
      } else {
        thevartab <- sapply(which.strata,
                            simplify=FALSE,
                            function(x) {
                              cbind(
                                  as.data.frame(theresults[which.vars,c(which.stats[!(which.stats=="p")]),x]),
                                  " " = format(Signif[which.vars,,x]))
                            })
      }
    }
    if (show.pvals & ("p" %in% dimnames(theresults)[["stat"]]) & !is.null(theresults)) {
      if (horizontal) {
        theftab <- ftabler(
            theresults[which.vars,which.stats,which.strata,drop=FALSE])
        thevartab <- theftab
      } else {
        thevartab <- sapply(
            which.strata,
            simplify=FALSE,
            function(x) {
              as.data.frame(theresults[which.vars,which.stats,x])
            })
      }
    }

    ##if(show.pvals&!("p"%in% dimnames(theresults)[["stat"]])& !is.null(theresults)) {
    ##  stop("You need to request p-values when calling xBalance.")
    ##} ##irrelevant now.

    if (!is.null(theoverall)) {
      nc <- length(theresults)/2
      latex.annotation <- NULL
      ## paste("\\\\ \\hline Overall",
      ##       paste("\\multicolumn{",nc,"}{c}{",preSig,"}"),
      ##       paste("\\multicolumn{",nc,"}{c}{",postSig,"}"),
      ##       sep=" & ")
      if (show.signif.stars) {
        ChiSignif <- signifier(theoverall[which.strata,"p.value"])

        theoveralltab <- cbind(format(theoverall[which.strata,],digits=DIGITS),format(ChiSignif))
        names(theoveralltab)[4]<-" "
      }
      theoveralltab<-format(theoverall[which.strata,],digits=DIGITS)
    } else {
      theoveralltab<-NULL
    }

    if (printme) {
      ## RItools:::print.ftable(thevartab,justify.labels="center",justify.data="right") ##doesn't seem to help the alignment problem
      if (!is.null(theresults)){ print(thevartab) }
      if (!is.null(theoverall) && print.overall) {
        cat("---Overall Test---\n")
        print(theoveralltab)
        if (show.signif.stars&!show.pvals) {
          if (!is.null(theresults)) {
            thelegend<-attr(Signif, "legend") ##if we are showing thevartab use the legend from that object
          }
          if (is.null(theresults) & !is.null(theoverall)) {
            thelegend<-attr(ChiSignif,"legend") ##use legend from the overall object if only showing that one
          }
          cat("---\nSignif. codes: ", thelegend, "\n")
        }
      }
      ##cat(paste("n = ", n, ", k = ", k,
      ##          "\nresidual sd = ", fround(summ$sigma, digits), ", R-Squared = ", fround(summ$r.squared, 2),
      ##          "\n", sep = ""))
    } else {
      list(vartable=thevartab,overalltable=theoveralltab)}
  }
              )
}
