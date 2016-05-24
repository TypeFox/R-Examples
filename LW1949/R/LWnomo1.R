#' Recreate Litchfield and Wilcoxon's Nomograph No. 1
#'
#' Recreate Litchfield and Wilcoxon's (1949) nomograph to estimate the
#' contribution to the chi-squared from the expected percent effect and
#' the observed  minus the expected percent effect.
#' @export
#' @import
#'   plotrix
#' @details
#' Use the nomograph by laying a straight edge from the expected percent effect
#' in the first scale to the observed (corrected, if necessary) minus the
#' expected percent effect in the second scale and reading the point where the
#' straight edge crosses the third scale as the contribution.
#'
#' The formula behind the nomograph is
#' (observed - expected)^2 / (100 * expected)
#' @param values
#'   A logical scalar indicating whether values should be output.
#' @return
#'   If \code{values} is TRUE, a list of length four, with the x and y
#'   coordinates and the corresponding values (all displayed in the log10
#'   scale) of the end points of the three scales.  Information is provided
#'   twice for the first scale, once for the left tick marks and once for the
#'   right tick marks.
#' @import
#'   graphics
#' @references
#'   Litchfield, JT Jr. and F Wilcoxon.  1949.
#'     A simplified method of evaluating dose-effect experiments.
#'     Journal of Pharmacology and Experimental Therapeutics 96(2):99-113.
#'     \href{http://jpet.aspetjournals.org/content/96/2/99.abstract}{[link]}.
#' @examples
#' LWnomo1()

LWnomo1 <- function(values=FALSE) {

  bigtix <- function(x, fudge=10, roundingto=c(1, 2, 5)) {
    onedigit <- signif(x, 1) - round(x, fudge) == 0
    gooddigit <- substring(format(signif(x, 1), sci=TRUE), 1, 1) %in% roundingto
    onedigit & gooddigit
    }

  # 1st scale, ep,
  #   expected % on log scale, log10(ep)
  ep2l <- c(
    seq(50,     80,    5),    seq(82,    90,   2),    seq(91,    95,    1),
    seq(95.5,   98,    0.5),  seq(98.2,  99,   0.2),  seq(99.1,  99.5,  0.1),
    seq(99.55,  99.8,  0.05), seq(99.82, 99.9, 0.02), seq(99.91, 99.95, 0.01),
    seq(99.955, 99.98, 0.005))
  ep1l <- rev(100-ep2l)
  ep1l. <- sort(unique(c(range(ep1l), ep1l[bigtix(ep1l)])))
  ep2l. <- rev(100 - ep1l.)

  # 3rd scale, chicont,
  #   100 times the contrib. to the chi-squared divided by n
  #   on the log scale, log10(100*contrib/n), where n is the total number
  chicontl <- 100*
c(seq(0.001, 0.002, 0.0002), seq(0.0025, 0.005, 0.0005),seq(0.006, 0.01, 0.001),
  seq(0.012, 0.02,  0.002),  seq(0.025,  0.05,  0.005), seq(0.06,  0.1,  0.01),
  seq(0.12,  0.2,   0.02),   seq(0.25,   0.5,   0.05),  seq(0.6,   1,    0.1),
  seq(1.2,   2,     0.2))

  chicontladj <- chicontl/100
  chicontl. <- sort(unique(c(range(chicontl), chicontl[bigtix(chicontl)])))
  chicontladj. <- chicontl./100

  # 2nd scale, opmep,
  #   observed minus expected % on log scale times 2, 2*log10|op - ep|
  # range of 2nd scale, as the sum of the ranges of the 1st and 3rd scales
  # ep + chicont = opmep
  opmeprange <- 10^((log10(c(0.02, 50)) + log10(c(0.1, 200)))/2)
  opmepladj <- sort(unique(c(opmeprange,
    seq(0.05, 0.1, 0.01), seq( 0.12, 0.2, 0.02), seq( 0.25, 0.5, 0.05),
    seq(0.6,  1,   0.1),  seq( 1.2,  2,   0.2),  seq( 2.5,  5,   0.5),
    seq(6,   10,   1),    seq(12,   20,   2),    seq(25,   50,   5),
    seq(60, 100,  10))))
  opmepl <- 2*log10(opmepladj)
  opmepladj. <- sort(unique(c(range(opmepladj), opmepladj[bigtix(opmepladj)])))
  opmepl. <- 2*log10(opmepladj.)

  par(xaxs="i", yaxs="i", mar=c(1, 1, 4, 1), las=1)
  plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
  axis(2, pos=0.1,  at=rescale(log10(ep1l),  0:1), labels=FALSE, tck=-0.01)
  axis(2, pos=0.1,  at=rescale(log10(ep1l.), 0:1), labels=round(rev(ep2l.), 2))
  axis(2, pos=0.152, at=rescale(log10(ep1l.), 0:1), labels=round(ep1l., 2),
    tick=FALSE, hadj=0)
  axis(2, pos=0.5, at=rescale(opmepl,  0:1)[-1], labels=FALSE, tck=-0.01)
  axis(2, pos=0.5, at=rescale(opmepl., 0:1)[-1],
    labels=round(opmepladj., 3)[-1])
  axis(2, pos=0.9, at=rescale(log10(chicontl),  0:1), labels=FALSE, tck=-0.01)
  axis(2, pos=0.9, at=rescale(log10(chicontl.), 0:1),
    labels=round(chicontladj., 4))
  mtext(c("Expected\n% effect", "Observed minus\nexpected % effect",
    "(Chi)\U00B2\nfor samples\nof one"), side=3, at=c(0.1, 0.5, 0.9), line=1)

  if(values) {
    scale1l <- data.frame(x= c(0.1, 0.1), y=0:1, values=c(99.98, 50))
    scale1r <- data.frame(x= c(0.1, 0.1), y=0:1, values=c(0.02, 50))
    scale2  <- data.frame(x= c(0.5, 0.5), y=0:1, values=c(0.045, 100))
    scale3  <- data.frame(x= c(0.9, 0.9), y=0:1, values=c(0.001, 2))
    out <- list(scale1l=scale1l, scale1r=scale1r, scale2=scale2, scale3=scale3)
    return(out)
  }
}
