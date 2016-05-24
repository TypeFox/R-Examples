#' Same sex marriage in the US
#'
#' Changes in the legality of same sex marriage in the United States over time.
#'
#' @source \url{http://www.nytimes.com/interactive/2015/03/04/us/gay-marriage-state-by-state.html}
#' @format A data frame with columns:
#' \describe{
#'  \item{State}{State Abbreviation}
#'  \item{Status}{Legal status. Either \code{bbs} meaning banned by statute,
#'    \code{nl} meaning not legal, \code{legal}, \code{bbca} meaning banned by
#'    constitutional ammendment, or \code{dis} meaning disputed.}
#'  \item{Year}{Year status went into effect.}
#' }
#' @examples
#' \dontrun{
#'  ssm
#' }
"ssm"

#' Postal Abbreviations for The United States of America
#'
#' @examples
#' \dontrun{
#'  usa_abb
#' }
"usa_abb"

#' Postal Abbreviations for Canada
#'
#' @examples
#' \dontrun{
#'  canada_abb
#' }
"canada_abb"

#' Postal Abbreviations for Mexico
#'
#' @examples
#' \dontrun{
#'  mexico_abb
#' }
"mexico_abb"

#' Production and farm value of maple products in Canada
#'
#' @source Statistics Canada. Table 001-0008 - Production and farm value of
#'  maple products, annual. \url{http://www5.statcan.gc.ca/cansim/}
#' @format A data frame with columns:
#' \describe{
#'  \item{Year}{A value between 1924 and 2015.}
#'  \item{Syrup}{Maple products expressed as syrup, total in thousands of gallons.}
#'  \item{CAD}{Gross value of maple products in thousands of Canadian dollars.}
#'  \item{Region}{Postal code abbreviation for territory or province.}
#' }
#' @examples
#' \dontrun{
#'  maple
#' }
"maple"

#' Monthly milk production in Canada
#'
#' @source Statistics Canada. Table 003-0011 - Milk production and utilization,
#' monthly. \url{http://www5.statcan.gc.ca/cansim/}
#' @format A data frame with columns:
#' \describe{
#'  \item{Year}{A value between 1976 and 2015.}
#'  \item{Month}{A value between 1 and 12.}
#'  \item{Region}{Postal code abbreviation for territory or province.}
#'  \item{Kiloliters}{Milk sold off farms in kiloliters.}
#' }
#' @examples
#' \dontrun{
#'  milk
#' }
"milk"
