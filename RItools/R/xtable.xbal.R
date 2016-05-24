##' This function uses the \code{\link[xtable]{xtable}} package
##' framework to display the results of a call to
##' \code{\link{xBalance}} in LaTeX format. At the moment, it ignores
##' the omnibus chi-squared test information.
##'
##' The resulting LaTeX will present one row for each variable in the
##' formula originally passed to \code{\link{xBalance}}, using the
##' variable name used in the original formula. If you wish to have
##' reader friendly labels instead of the original variables names,
##' see the code examples below.
##'
##' To get decimal aligned columns, specify \code{align=c("l",
##' rep(".", <ncols>))}, where \code{<ncols>} is the number of columns
##' to be printed, in your call to \code{xtable}.  Then use the
##' \code{dcolumn} package and define \samp{'.'} within LaTeX: add the
##' lines \code{\\usepackage\{dcolumn\}} and
##' \code{\\newcolumntype\{.\}\{D\{.\}\{.\}\{2.2\}\}} to your LaTeX
##' document's preamble.
##' @title An \code{xtable} method for \code{xbal} objects
##' @param x An object resulting from a call to
##'   \code{\link{xBalance}}.
##' @param caption See \code{\link[xtable]{xtable}}.
##' @param label See \code{\link[xtable]{xtable}}.
##' @param align See \code{\link[xtable]{xtable}}. Our default (as of
##'   version 0.1-7) is right-aligned columns; for decimal aligned
##'   columns, see details, below.
##' @param digits See \code{\link[xtable]{xtable}}. Default is 2.
##' @param display See \code{\link[xtable]{xtable}}.
##' @param col.labels Labels for the columns (the test
##'   statistics). Default are come from the call to
##'   \code{\link{print.xbal}}.
##' @param ... Other arguments to \code{\link{print.xbal}}.
##' @return This function produces an \code{xtable} object which can
##'   then be printed with the appropriate \code{print} method (see
##'   \code{\link[xtable]{print.xtable}}).
##' @export
##' @import xtable
##' @examples
##' data(nuclearplants)
##' require(xtable)
##'
##' # Test balance on a variety of variables, with the 'pr' factor
##' # indicating which sites are control and treatment units, with
##' # stratification by the 'pt' factor to group similar sites
##' xb1 <- xBalance(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'                 strata = data.frame(unstrat = factor(character(32)),
##'                 pt = factor(nuclearplants$pt)),
##'                 data = nuclearplants,
##'                 report = c('adj.means', 'adj.mean.diffs',
##'                            'std.diffs', 'z.scores',
##'                            'chisquare.test', 'p.values'))
##'
##' xb1.xtab <- xtable(xb1) # This table has right aligned columns
##'
##' # Add user friendly names in the final table
##' rownames(xb1.xtab) <- c("Date", "Application to Contruction Time",
##' "License to Construction Time", "Net Capacity", "Northeast Region", "Cooling Tower",
##' "Babcock-Wilcox Steam", "Cumlative Plants")
##'
##' print(xb1.xtab,
##'       add.to.row = attr(xb1.xtab, "latex.add.to.row"),
##'       hline.after = c(0, nrow(xb1.xtab)),
##'       sanitize.text.function = function(x){x},
##'       floating = TRUE,
##'       floating.environment = "sidewaystable")
xtable.xbal <- function(x,
                        caption = NULL,
                        label = NULL,
                        align = c("l",rep("r",ncol(xvardf))),
                        digits = 2,
                        display = NULL,
                        col.labels=NULL,
                        ...) {
  ##By default use decimal alignment, which will require the dcolumn package in latex and an appropriate column definition like:
  ##\newcolumntype{.}{D{.}{.}{2.2}}
  ##Here is an example which works
  ##xb1<-xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
  ##         strata=data.frame(unstrat=factor(character(32)),
  ##           pt=factor(nuclearplants$pt)),
  ##         data=nuclearplants,
  ##         report=c("adj.means","adj.mean.diffs",'std.diffs', 'z.scores', 'chisquare.test','p.values'))
  ##
  ##junk<-xtable(xb1)
  ##print(junk,add.to.row=attr(junk,"latex.add.to.row"),hline.after=c(0,nrow(junk)),sanitize.text.function=function(x){x},floating=TRUE,floating.environment="sidewaystable")

  xprint <- flatten.xbalresult(x)
  numstrata<-dim(x$results)[3]
  latex.annotation <- attr(xprint, "latex.annotation")
  xvardf<-xprint$vartable

  if (!is.null(col.labels))
    names(xvardf) <- col.labels


  ##call xtable on the resulting data.frame
  vartab <- xtable(xvardf,caption=caption, label=label, digits=digits,align=align,display=display,col.labels=col.labels,...) ##NextMethod("xtable",xvardf)
  structure(vartab,
            latex.add.to.row=list(pos=list(-1),command=latex.annotation),
            hline.after=c(0,nrow(xvardf)))

}
