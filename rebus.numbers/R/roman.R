#' Roman numerals
#'
#' Match roman numerals.
#' @param lo A non-negative integer. Minimum number of repeats, when grouped.
#' @param hi positive integer. Maximum number of repeats, when grouped.
#' @examples
#' # Constant form and character class
#' ROMAN
#' roman()
#'
#' x <- c("MMMDCCCXLVIII", "MMMCMDCCCXLVIIV")
#' rx <- rebus.base::exactly(roman())
#' grepl(rx, x)
#' @name roman
#' @seealso \code{\link[utils]{as.roman}}
#' @importFrom utils as.roman
#' @export
ROMAN <- case_insensitive(
    group(
    repeated("M", 0, 3) %R%
    optional(or1(as.character(as.roman(seq.int(100, 900, 100))))) %R%
    optional(or1(as.character(as.roman(seq.int(10, 90, 10))))) %R%
    optional(or1(as.character(as.roman(1:9))))
  )
)


#' @rdname roman
#' @export
roman <- function(lo, hi)
{
  repeated(ROMAN, lo, hi, char_class = FALSE)
}
