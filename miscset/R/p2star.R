#' @name p2star
#' @author Sven E. Templer
#' @title P Value Significance Level Indicator
#' @description 
#' Transform p-values to character (e.g. stars) indicators
#' by significance levels with the function \link{symnum}.
#' @param p Vector with p values
#' @param breaks The breaks from min (0) to max (1).
#' @param symbols Symbols to use for values between breaks from min to max.
#' @examples
#' #
#' 
#' p2star(c(1e-5,.1,.9))
#' 
#' #

#' @rdname p2star
#' @export p2star
p2star <- function(
  p, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
  symbols = c("***", "**", "*", ".", "n.s.")) 
{
  as.character(symnum(p, breaks, symbols))
  #cut(p, breaks, labels, include.lowest = T)
}

# benchmarks:
# symnum: .84 (length(p)=1e3) 1.6 (1e4) 11.6 (1e5) 201.7 (1e6) ms
# cut   : .42 (length(p)=1e3) 3.0 (1e4) 30.9 (1e5) 307.1 (1e6) ms
