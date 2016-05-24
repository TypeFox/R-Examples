#' An S4 class to represent the top-level option model
#' 
#' The \code{ecop} class serves as an object-oriented container for
#' the option pricing model. It does have a specific purpose at the moment -
#' that is, to produce all the data for the charts of the paper, based on
#' CBOE data structure. Therefore, user may not find it general enough.
#' That probably will be the case for the time being until more popularity calls
#' for a more generic container.
#'
#' @slot call      the match.call slot
#' @slot conf      list, configuration
#' @slot key       character
#' @slot symbol    character
#' @slot datadate  Date
#' @slot days      numeric, days between datadate and expiry date
#' @slot ttm       numeric, time to maturity in days/365
#' @slot int_rate  numeric
#' @slot div_yield numeric
#' @slot put_data  the put data of ecop.opt class
#' @slot call_data the call data of ecop.opt class
#' @slot put_conf  list, the put configuration
#' @slot call_conf list, the call configuration
#'
#' @include ecop-opt-class.R
#' @keywords ecop
#' 
#' @author Stephen H-T. Lihn
#'
#' @exportClass ecop
setClass("ecop",
         representation(call = "call",
                        conf = "list",
                        key = "character",
                        symbol = "character",
                        datadate = "Date",
                        days = "numeric",
                        ttm = "numeric",
                        int_rate = "numeric",
                        div_yield = "numeric",
                        put_data = "ecop.opt",
                        call_data = "ecop.opt",
                        put_conf = "list",
                        call_conf = "list"
                       ),
          prototype(call = call("ecop"),
                    conf = list(),
                    key = "",
                    symbol = "",
                    datadate = as.Date("1970-01-01"),
                    days = NaN,
                    ttm = NaN,
                    int_rate = NaN,
                    div_yield = NaN,
                    put_data = NULL,
                    call_data = NULL,
                    put_conf = list(),
                    call_conf = list()
                   )
)
### <---------------------------------------------------------------------->
