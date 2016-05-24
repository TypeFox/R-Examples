#' setClass for ecdq class
#' 
#' setClass for ecdq class, the quantile generator
#'
#' @slot call the match.call slot
#' @slot xseg.from,xseg.to numeric vectors. The from and to for each x segment.
#' @slot cseg.from,cseg.to numeric vectors. The from and to for each cdf segment.
#' @slot cseg.min,cseg.max numeric. The min and max of cdf segments.
#' @slot N_seg numeric. Number of segments.
#' @slot cdf.fit A vector of \code{lm} object, one for each segment.
#' @slot x_left_tail,x_right_tail numeric. The starting x of left and right tails.
#' @slot fit.left,fit.right objects of \code{lm} class for fitting the tails.
#' @slot conf list of miscelaneous configurations. For debugging purpose.

#' @keywords ecdq class constructor
#' @include ecd-numericMpfr-class.R
#'
#' @exportClass ecdq
setClass("ecdq",
         representation(call = "call",
                        xseg.from = "numericMpfr",
                        xseg.to   = "numericMpfr",
                        cseg.from = "numericMpfr",
                        cseg.to   = "numericMpfr",
                        cseg.min  = "numericMpfr",
                        cseg.max  = "numericMpfr",
                        N_seg     = "numericMpfr",
                        cdf.fit = "list",
                        x_left_tail  = "numericMpfr",
                        x_right_tail = "numericMpfr",
                        fit.left  = "lm",
                        fit.right = "lm",
                        conf = "list"),
          prototype(call = call("ecdq"),
                    xseg.from = NaN,
                    xseg.to   = NaN,
                    cseg.from = NaN,
                    cseg.to   = NaN,
                    cseg.min  = NaN,
                    cseg.max  = NaN,
                    N_seg     = NaN,
                    cdf.fit = list(),
                    x_left_tail  = NaN,
                    x_right_tail = NaN,
                    fit.left  = "lm",
                    fit.right = "lm",
                    conf = list())
)
### <---------------------------------------------------------------------->
