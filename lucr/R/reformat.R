#'@title Format numeric or integer values as currency values
#'
#'@description a formatter that lets you take numeric or integer values (12000) and convert them
#'to strings that are currency-formatted ($12,000). Full control is available over the currency symbol,
#'the size of delimited groups, the sign used \emph{to} delimit groups, and decimal placement.
#'
#'@param x a numeric or integer vector containing values you want to currency-ify
#'
#'@param currency_symbol the symbol that identifies the currency. "£" by default.
#'
#'@param symbol_first whether the symbol goes at the beginning (TRUE) or end (FALSE) of the generated
#'value. TRUE by default.
#'
#'@param group_size the size of delimited groups (2, 3, or 4 digits, say). Set to 3 by default.
#'
#'@param group_delim the delimiter for each group.
#'
#'@param decimal_size the number of digits after the decimal place. 2 by default but can be more (the
#'Japanese Yen, for example, can go down to one \emph{rin}, which is a thousandth of a Yen).
#'
#'@param decimal_delim the delimiter to use for sub-unit, decimal values. A period by default.
#'
#'@examples
#'to_currency(120000.03)
#'#[1] "£120,000.03"
#'
#'@export
to_currency <- function(x, currency_symbol = "\u00A3", symbol_first = TRUE, group_size = 3,
                        group_delim = ",", decimal_size = 2, decimal_delim = "."){
  
  if(is.integer(x)){
    return(currency_format_(x = format(x, scientific = 300), currency_symbol, symbol_first,
                            group_size, decimal_delim, group_delim))
  } else if(is.numeric(x)){
    return(currency_format_(x = format(x, nsmall = decimal_size, scientific = 300), currency_symbol, symbol_first,
                            group_size, decimal_delim, group_delim))
  }
  stop("x must be a numeric or integer vector")
}

#'@title Convert currency-formatted strings into numeric values
#'@description takes a vector of strings formatted as amounts of money ("$12,329.34")
#'and reformats them as numeric values (12329.34).
#'
#'@param x a vector of strings, formatted as money amounts. See \code{\link{to_currency}}.
#'
#'@param decimal_delim the character used to delimit the decimal amount. Set to "." by default.
#'
#'@examples
#'from_currency("£1,249.34")
#'# [1] 1249.34
#'@export
from_currency <- function(x, decimal_delim = "."){
  currency_unformat_(x, decimal_delim)
}
