#' Ruby-like string interpolation
#'
#' To interpolate any R object into a string, it can be wrapped in \code{#\{...\}}
#' and within a string. This is equivalent to evaluating whatever is in this block
#' within the current environment (see \code{base::environment} when called with
#' no arguments).
#'
#' @param ... arbitrary many character vectors, that will be concatenated
#'   using the \code{sep} parameter as delimiter. Any non-length 1 vectors
#'   will also be collapsed using the \code{collapse} parameter.
#' @param envir environment. Where to evaluate any interpolated strings,
#'   by default \code{parent.frame()}, the current local environment.
#' @param sep character. The delimiter to use to concatenate passed
#'   strings, by default the empty string \code{''}.
#' @param collapse character. The delimiter to use to concatenate passed
#'   character vectors of length greater than 1, by default the empty string
#'  \code{''}.
#' @export
#' @examples
#' x <- 1; y <- 2
#' cat(pp("1 + #{x} = #{1 + x}")) # "1 + 1 = 2"
#' cat(pp(1:3, " are numbers")) # "123 are numbers"
#' cat(pp(1:3, " are numbers", collapse = " and "))
#' # "1 and 2 and 3 are numbers"
#' cat(pp("x = #{x}", "y = #{y}", "x + y = #{x + y}", sep = ", "))
#' # x = 1, y = 2, x + y = 3
pp <- function(..., envir = parent.frame(), sep = '', collapse = '') {
  string <- list(...)
  if (length(string) > 1)
    return(paste(sapply(string,
      function(s) { pp(s, envir = envir, sep = sep, collapse = collapse) }
    ), collapse = sep))
  string <- string[[1]]
  if (length(string) > 1)
    return(paste(sapply(string,
      function(s) { pp(s, envir = envir, sep = sep, collapse = collapse) }
    ), collapse = collapse))
  regex <- gregexpr('#\\{([^\\}]+)\\}', string, perl=TRUE)
  starts <- attr(regex[[1]], 'capture.start')
  lengths <- attr(regex[[1]], 'capture.length')
  buildstr <- ''
  last <- 1
  for(i in 1:length(attr(regex[[1]], 'capture.start'))) {
    buildstr <- append(buildstr,
      c(substr(string, last, starts[i] - 3),
        eval(parse(text = substr(string, starts[i], starts[i] + lengths[i] - 1)),
        envir = envir)
       ))
    
    last <- starts[i] + lengths[i] + 1
  }
  buildstr <- append(buildstr, substr(string, last, nchar(string)))
  paste(buildstr, collapse = '')
}
