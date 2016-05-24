#' Commodity Futures Data from 2003 to 2005.
#' 
#' This data frame contains the futures backtest results of 28 commodities from 
#' 2003-01-01 to 2005-12-30. The data frame contains daily positions and profits
#' for each commodity.
#' 
#' The commodity data frame contains three layers of complexity. The first layer
#' is the strategy layer. The \code{strategy} column contains three different 
#' strategies that divide the whole data frame into three part. Each strategy is
#' further divided into different substrategies, and all these substrategies are
#' contained in the \code{substrategy} column. For each substrategy, different 
#' portfolios are formed regularly and overlapping with each other, so each 
#' substrategy is divided into different overlapping portfolios, contained in
#' the \code{portfolio} column.
#' 
#' @format A data frame with 11 variables 
#' \itemize{ 
#'   \item name = The name of each commodity.
#'   \item id = The specific ID for each commodity.
#'   \item date = The individual trading date.
#'   \item sector = The sector a commodity.
#'   \item portfolio = The portfolio number the holding belongs to.
#'   \item strategy = The strategy under which substrategies and portfolios are
#'         established. The same commodity can belong to different strategies at
#'         the same time.
#'   \item substrategy = The substrategy under which portfolios are established.
#'         Substrategy is only a finer division of strategies.
#'   \item gmv = The gross market value of the commodity on the trading date.
#'   \item nmv = The net market value of the commodity on the trading date.
#'   \item pnl = The adjusted P&L of a commodity on the trading date.
#'   \item contract = The number of contracts of a commodity on that trading 
#'         date.}
#' @docType data
#' @name commodity
#' @keywords commodity backtest data
NULL