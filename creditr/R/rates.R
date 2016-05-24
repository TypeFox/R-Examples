#' LIBOR rates from 2004-01-01 to 2015-08-03
#' 
#' This data frame is created by \code{build_rates} in \code{creditr} package to
#' calculate the CDS pricing. It covers three currencies: USD, EUR, and JPY. The
#' interest rates date from 2004-01-01 to 2014-08-23. Rates on holidays and 
#' weekends are available as well as business days.
#' 
#' @format A data frame with 194378 observations on the following 4 variables. 
#'   \itemize{ \item date = a date (Date object) \item currency = a character 
#'   containing \code{USD} \code{EUR} \code{JPY} \item expiry = a character 
#'   containing \code{1M} \code{2M} \code{3M} \code{6M} \code{9M} \code{1Y} 
#'   \code{2Y} \code{3Y} \code{4Y} \code{5Y} \code{6Y} \code{7Y} \code{8Y} 
#'   \code{9Y} \code{10Y} \code{12Y} \code{15Y} \code{20Y} \code{30Y} \item rate
#'   = a numeric vector. The LIBOR rate. }
#' @details The source of the interest rates in \code{rates.RData} is from 
#'   \url{https://www.markit.com/} and 
#'   \url{http://research.stlouisfed.org/fred2/}. When a user is calculating CDS
#'   using the \code{CDS} package, the package calls \code{get_rates} to get the
#'   needed interest rates; \code{get_rates} then calls the \code{rates.RData} 
#'   for these interest rates. If a date is unavailable in \code{rates.RData}, 
#'   then the package calls other functions to get the needed interest rates 
#'   from the internet. The \code{rates.RData} is created and stored in the 
#'   package for the users' convenience: getting interest rates from the 
#'   Internet may fail due to the internet connection problem and may be very 
#'   slow. Also, the user can build its own updated local \code{rates.RData} by 
#'   using \code{build_rates}. For more explanation on the usage of 
#'   \code{build_rates}, please see \bold{See Also}.
#'   
#'   Also, please notice that in the \code{rates.RData}, the \code{rate} is not 
#'   the interest rate on the \code{date} listed in the same row as the 
#'   \code{rate}. For example, in the first row, the \code{rate} is 0.001550; 
#'   the \code{date} in that row is 2014-08-21. But this does not mean that on 
#'   2014-08-21, the interest rate is 0.001550; instead, 0.001550 is the 
#'   interest rate on 2014-08-20. This regulation of using the interest rate on 
#'   the previous business day of CDS trading date is set by ISDA Standard 
#'   Model. The Model says that, if a trader buy a CDS on 2014-08-21, then when 
#'   calculating the pricing of the CDS, she should use the interest rate of 
#'   2014-08-20, which is 0.001550. This may be confusing for people who are yet
#'   unfamiliar with CDS. We design \code{rates.RData} in this way because it is
#'   easier for other functions in the \code{CDS} package to use this data frame
#'   for calculation.
#'   
#'   \code{rates.RData} covers holidays, weekend and business days. As is set by
#'   ISDA Standard Model and introduced above, we use the previous business 
#'   day's interest rate for CDS pricing on a certain trading date. Therefore, 
#'   if the user is buying a CDS on Saturday, she should use the interest rate 
#'   of last Friday; if she is buying a CDS on Sunday or Monday, she should 
#'   still use the interest rate of last Friday, because last Friday is the 
#'   previous business day of the CDS trading date. When it comes to holidays, 
#'   we still choose the previous business day for interest rate. For example, 
#'   if a trader is buying a CDS on 2014-07-05, then she should use the interest
#'   rate of 2014-07-03, because 2014-07-04 is a national holiday and 2014-07-03
#'   is the previous business day of trading date.
#'   
#'   Also, please notice that in \code{rates.RData}, a currency's type of 
#'   expiries generally stays the same along the time, but not always. For 
#'   example, we check the expiry type of USD in \code{rates.RData}: for 
#'   2004-01-01, there are two types of expirty: 1M, 3Y; for 2006-01-01, there 
#'   are five types of expiry: 1M, 2M, 3M, 6M, 1Y; for 2007-01-01, there are 18 
#'   types of expiry; for 2014-01-01, there are 19 types of expiry. Therefore, 
#'   for a trader, she should always keep in mind to update her knowledge of the
#'   expiry types of her trading currency. For different currencies, the types 
#'   of expiries are often different.
#'   
#'   Finally, please notice that some of the data are missing for a certain 
#'   expiry of a currency in a short time. For example, some dates in 
#'   \code{rates.RData} do not have a \code{expiry} of 3Y for \code{USD}. This 
#'   is not likely to be caused by data error, since all these data are got from
#'   Markit and FRED. Users, however, should be aware of that some data seem 
#'   "missing".
#'   
#' @source \url{https://www.markit.com/} 
#'   \url{http://research.stlouisfed.org/fred2/}
#'   
#' @examples 
#' data(rates)
#' 
#' ## for JPY rates: 
#' rates[rates$currency == "JPY",]
#' 
#' ## for rates on a specific date, of a specific currency:
#' rates[rates$currency == "USD" & rates$date == "2005-10-01",]
#' 
#' @docType data
#' @name rates
#' @keywords datasets, interest rates
#' @seealso \code{\link{download_FRED}} \code{\link{download_markit}}
#'   \code{\link{build_rates}}
NULL