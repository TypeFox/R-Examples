##' Export R objects to several markup languages
##'
##' Convert an R object to an \code{ascii} object, which can then be
##' printed with asciidoc, txt2tags, reStructuredText, org, textile or
##' pandoc syntax.
##'
##' The nature of the generated output depends on the class of
##' \code{x}.  For example, \code{summary.table} objects produce a
##' bulleted list while \code{data.frame} objects produce a table of
##' the entire data.frame.
##'
##' Sometimes, arguments are not active, depending of the features
##' implemented in the markup language generated. All arguments are
##' active when asciidoc syntax is produced.
##'
##' The available method functions for \code{ascii} are given by
##' \code{methods(ascii)}.  Users can extend the list of available
##' classes by writing methods for the generic function \code{ascii}.
##' All method functions should return an object of class
##' \code{"ascii"}.
##'
##' @aliases ascii package-ascii ascii.anova ascii.aov ascii.aovlist ascii.cast_df ascii.character ascii.coxph ascii.CrossTable ascii.data.frame ascii.default ascii.density ascii.describe ascii.describe.single ascii.factor ascii.freqtable ascii.ftable ascii.glm ascii.htest ascii.integer ascii.list ascii.lm ascii.matrix ascii.meanscomp ascii.mtable ascii.numeric ascii.packageDescription ascii.prcomp ascii.sessionInfo ascii.simple.list ascii.smooth.spline ascii.stat.table ascii.summary.aov ascii.summary.aovlist ascii.summary.glm ascii.summary.lm ascii.summary.prcomp ascii.summary.survfit ascii.summary.table ascii.survdiff ascii.survfit ascii.table ascii.ts ascii.zoo
##' @return This function returns an object of class
##' \code{"asciiTable"}, \code{"asciiList"} or \code{"asciiMixed"}.
##' @author David Hajage \email{dhajage@@gmail.com}
##' @keywords print
##' @rdname ascii
##' @export
##' @examples
##' data(esoph)
##' ascii(esoph[1:10,])
##' tab <- table(esoph$agegp, esoph$alcgp)
##' ascii(tab)
##' print(ascii(tab), type = "t2t")
##' print(ascii(tab), type = "rest")
##' print(ascii(tab), type = "org")
##' ascii(summary(tab))
##'
ascii <- function (x, ...) {
  UseMethod("ascii")
}
