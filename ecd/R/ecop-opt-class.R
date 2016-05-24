#' An S4 class to represent the option data and model calculation
#' 
#' The \code{ecop.opt} class serves as an object-oriented container for
#' the type-specific (p or c) option data.
#'
#' @slot call the match.call slot
#' @slot otype character, option type
#' @slot range.from numeric, starting price range
#' @slot range.to numeric, ending price range
#' @slot momentum numeric, momentum for tranlation (T) operator
#' @slot epsilon numeric, asymptotic premium
#' @slot k_cusp numeric, the suggested cusp location for poly fit of prices
#' @slot ecldOrEcd the ecld/ecd class to calculate theoretical values in local regime
#' @slot S underyling price, this can be overriden by conf
#' @slot S_raw underyling price (before override)
#' @slot strike strike price
#' @slot k log-strike price
#' @slot V_last last option price
#' @slot V_bid bid option price
#' @slot V_ask ask option price
#' @slot V finalized option price (likely mid-point)
#'
#' @include ecd-ecldOrEcd-class.R
#' @keywords ecop
#' 
#' @author Stephen H. Lihn
#'
#' @exportClass ecop.opt
setClass("ecop.opt",
         representation(call = "call",
                        otype = "character",
                        range.from = "numeric",
                        range.to = "numeric",
                        momentum = "numeric",
                        epsilon = "numeric",
                        k_cusp = "numeric",
                        ecldOrEcd = "ecldOrEcd",
                        S = "numeric",
                        S_raw = "numeric",
                        strike = "numeric",
                        k = "numeric",
                        V_last = "numeric",
                        V_bid = "numeric",
                        V_ask = "numeric",
                        V = "numeric"
                    ),
          prototype(call = call("ecop"),
                    otype = "",
                    range.from = NaN,
                    range.to = NaN,
                    momentum = NaN,
                    epsilon = NaN,
                    k_cusp = NaN,
                    ecldOrEcd = NULL,
                    S = NaN,
                    S_raw = NaN,
                    strike = NaN,
                    k = NaN,
                    V_last = NaN,
                    V_bid = NaN,
                    V_ask = NaN,
                    V = NaN)
)
### <---------------------------------------------------------------------->
