#' @details \code{break_even} contains detailed run times for the analysis of
#' break-even points between the various implementations.
#'
#' @format A data frame with 5 columns:
#'
#' \describe{
#'   \item{\code{prob}}{
#'     A description of the probability distribution used.
#'     See \code{data_raw/benchmark.R} for details.
#'   }
#'   \item{\code{expr}}{
#'     Function name without the \code{sample_int_} prefix.
#'   }
#'   \item{\code{time}}{
#'     Run time in nanoseconds, as measured by
#'     \code{\link[microbenchmark]{microbenchmark}}
#'   }
#'   \item{\code{r}}{
#'     Ratio between the \code{size} and \code{n} arguments.
#'   }
#'   \item{\code{n}}{
#'     The \code{n} argument.
#'   }
#' }
#'
#' @rdname timings
#' @include timings.R
#' @name break_even
#' @export
#' @examples
#' head(break_even)
NULL
