#' rebus: Regular Expression Builder, Um, Something
#' 
#' Build regular expressions in a human readable way.
#' 
#' Regular expressions are a very powerful tool, but the syntax is terse enough 
#' to be difficult to read. This makes bugs easy to introduce, and hard to 
#' find. This package contains functions to make building regular expressions
#' easier.
#' @docType package
#' @name rebus
#' @aliases rebus rebus-package
#' @seealso \code{\link[base]{regex}} and \code{\link[base]{regexpr}}
#' The `stringr` and `stringi` packages provide tools for matching 
#' regular expressions and nicely complement this package.
#' \url{http://www.regular-expressions.info} has good advice on using
#' regular expression in R.  In particular, see
#' \url{http://www.regular-expressions.info/rlanguage.html} and
#' \url{http://www.regular-expressions.info/examples.html}
#' \url{https://www.debuggex.com} is a visual regex debugging and testing site.
#' @examples
#' ### Match a hex colour, like `"#99af01"`
#' # This reads *Match a hash, followed by six hexadecimal values.*
#'   
#' "#" %R% hex_digit(6)    
#' 
#' # To match only a hex colour and nothing else, you can add anchors to the 
#' # start and end of the expression.
#' 
#' START %R% "#" %R% hex_digit(6) %R% END
#' 
#' ### Simple email address matching. 
#' # This reads *Match one or more letters, numbers, dots, underscores, percents, 
#' # plusses or hyphens. Then match an 'at' symbol. Then match one or more letters, 
#' # numbers, dots, or hyphens. Then match a dot. Then match two to four letters.*
#'   
#' one_or_more(char_class(ASCII_ALNUM %R% "._%+-")) %R%
#'   "@@" %R%
#'   one_or_more(char_class(ASCII_ALNUM %R% ".-")) %R%
#'   DOT %R%
#'   ascii_alpha(2, 4)
#' 
#' ### IP address matching. 
#' # First we need an expression to match numbers between 0 and 255.  Both the 
#' # following syntaxes read *Match two then five then a number between zero and 
#' # five.  Or match two then a number between zero and four then a digit. Or match 
#' # an optional zero or one followed by an optional digit folowed by a compulsory 
#' # digit.  Make this a single token, but don't capture it.*
#' 
#' # Using the %|% operator
#' ip_element <- group(
#'   "25" %R% char_range(0, 5) %|%
#'   "2" %R% char_range(0, 4) %R% ascii_digit() %|%
#'   optional(char_class("01")) %R% optional(ascii_digit()) %R% ascii_digit()
#' )
#' 
#' # The same again, this time using the or function
#' ip_element <- or(
#'   "25" %R% char_range(0, 5),
#'   "2" %R% char_range(0, 4) %R% ascii_digit(),
#'   optional(char_class("01")) %R% optional(ascii_digit()) %R% ascii_digit()
#' )
#' 
#' # It's easier to write using number_range, though it isn't quite as optimal 
#' # as handcrafted regexes.
#' number_range(0, 255, allow_leading_zeroes = TRUE)
#' 
#' # Now an IP address consists of 4 of these numbers separated by dots. This 
#' # reads *Match a word boundary. Then create a token from an `ip_element` 
#' # followed by a dot, and repeat it three times.  Then match another `ip_element`
#' # followed by a word boundary.*
#' 
#' BOUNDARY %R% 
#'   repeated(group(ip_element %R% DOT), 3) %R% 
#'   ip_element %R%
#'   BOUNDARY    
#' @author Richard Cotton \email{richierocks@@gmail.com}
#' @include imports.R
NULL
