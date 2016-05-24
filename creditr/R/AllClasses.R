#' CDS Class
#' 
#' Class definition for the \code{CDS-Class}
#' 
#' @slot name is the name of the reference entity. Optional.
#' @slot contract is the contract type, default SNAC
#' @slot RED alphanumeric code assigned to the reference entity. Optional.
#'   
#' @slot date is when the trade is executed, denoted as T. Default is
#'   \code{Sys.Date}.
#' @slot spread CDS par spread in bps.
#' @slot maturity date of the CDS contract.
#' @slot tenor of contract in number of years - 5, 3
#' @slot coupon quoted in bps. It specifies the payment amount from
#' @slot recovery in decimal. Default is 0.4.
#' @slot currency in which CDS is denominated.
#' @slot principal is the dirty \code{upfront} less the \code{accrual}.
#' @slot accrual is the accrued interest payment.
#' @slot pd is the approximate the default probability at time t given the
#'   \code{spread}.
#' @slot price
#' @slot upfront is quoted in the currency amount. Since a standard contract is
#'   traded with fixed coupons, upfront payment is introduced to reconcile the
#'   difference in contract value due to the difference between the fixed coupon
#'   and the conventional par spread. There are two types of upfront, dirty and
#'   clean. Dirty upfront, a.k.a. Cash Settlement Amount, refers to the market
#'   value of a CDS contract. Clean upfront is dirty upfront less any accrued 
#'   interest payment, and is also called the Principal.
#' @slot spread.DV01 measures the sensitivity of a CDS contract mark-to-market
#'   to a parallel shift in the term structure of the par spread.
#' @slot IR.DV01 is the change in value of a CDS contract for a 1 bp parallel
#'   increase in the interest rate curve. \code{IRDV01} is, typically, a much
#'   smaller dollar value than \code{spreadDV01} because moves in overall
#'   interest rates have a much smaller effect on the value of a CDS contract
#'   than does a move in the CDS spread itself.
#' @slot rec.risk.01 is the dollar value change in market value if the recovery
#'   rate used in the CDS valuation were increased by 1\%.
#'   
#' @name CDS, CDS-class
#' @aliases CDS, CDS-class
#' @docType class
#' @rdname CDS-class
#' @export

setClass("CDS",
         representation = representation(
           
           ## name stuff
           
           name      = "character",
           contract  = "character",
           RED       = "character",
           
           ## basic info
           
           date      = "Date",
           spread    = "numeric",
           maturity  = "Date",
           tenor     = "numeric",
           coupon    = "numeric",
           recovery  = "numeric",
           currency  = "character",
           notional  = "numeric",
           principal = "numeric",
           accrual   = "numeric",
           pd        = "numeric",
           price     = "numeric",
           
           ## calculated amount
           
           upfront     = "numeric",
           spread.DV01 = "numeric",
           IR.DV01     = "numeric",
           rec.risk.01 = "numeric" ),
         
         prototype = prototype(
           
           ## name stuff
           
           name     = character(),
           contract = character(),
           RED      = character(),
           
           ## basic info
           
           date      = character(),
           spread    = numeric(),
           maturity  = character(),
           tenor     = numeric(),
           coupon    = numeric(),
           recovery  = numeric(),
           currency  = character(),
           notional  = numeric(),
           principal = numeric(),
           accrual   = numeric(),
           pd        = numeric(),
           price     = numeric(),
           
           ## calculated amount
           
           upfront     = numeric(),
           spread.DV01 = numeric(),
           IR.DV01     = numeric(),
           rec.risk.01 = numeric()
         )
)