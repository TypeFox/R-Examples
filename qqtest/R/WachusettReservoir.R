#' Storage, in millions of gallons daily per square mile of net land area, at the Wachusett Reservoir in Massacusetts - storage computed for each of several rates of draft (draft being a determined maintainable flow in 1,000s of gallons per square mile daily).
#'
#' First extensive published use of normal qqplots. Hazen uses a=1/2 to make the p values for the plots.  Hazen doesn'tplot zeros but has them contribute to the sample size. The context of use is in a study of the relation between the water storage provided in a reservoir on any stream and the quantity of water that can be continuously supplied by it.  To quote the paper: ... treat all the remaining variations on the basis of probabilities, using all data from a number of streams; and to study them in comparison with the normal law of error."
#' 
#' \code{qqtest(WachusettReservoir$draft800,dist="uniform", a=1/2,type="o")} will effect Hazen's original plot for a draft of 800,000 gallons per square mile daily.
#' 
#' \code{qqtest(WachusettReservoir$draft800,dist="normal", a=1/2, type="o")} will effect Hazen's normal qq plot for a draft of 800,000 gallons per square mile daily.
#'
#'
#' @format A data frame with 15 rows and 6 variates:
#' \describe{
#'   \item{draft100}{Computed storage, in millions of gallons per square mile of land area, given a draft of 100,000 gallons per square mile daily.}
#'   \item{draft200}{Computed storage, in millions of gallons per square mile of land area, given a draft of 200,000 gallons per square mile daily.}
#'   \item{draft400}{Computed storage, in millions of gallons per square mile of land area, given a draft of 400,000 gallons per square mile daily.}
#'   \item{draft600}{Computed storage, in millions of gallons per square mile of land area, given a draft of 600,000 gallons per square mile daily.}
#'   \item{draft800}{Computed storage, in millions of gallons per square mile of land area, given a draft of 800,000 gallons per square mile daily.}
#'   \item{draft1000}{Computed storage, in millions of gallons per square mile of land area, given a draft of 1,000,000 gallons per square mile daily.}
#' }
#' @source 
#' "Storage to be provided in impounding reservoirs for municipal water supply (with discussion)", 
#' Allen Hazen, Transactions of the American Society of Civil Engineers, Vol. 77, (1914), pp. 1539-1669.  
#'  
"WachusettReservoir"

WachusettReservoir <- data.frame(
						 draft100=c(0,0,0,0,0,
									  0,0,0,0,0,
									  0.2, 0.2,0.1,1.0,1.2),
						 draft200=c(0,0,0,2,0,
									 0,0,0,0,0,
									 3, 6, 3, 9, 4),
						 draft400=c(5, 2, 15, 23, 2,
									 11, 0.5, 3, 3, 3,
									 11, 27, 25, 34, 23),
						 draft600=c(17, 8, 52, 48, 14,
									 35, 10, 31, 20, 11,
									 23, 68, 57, 70, 51),
						 draft800=c(29, 14, 95, 78, 45,
									 59, 35, 68, 53, 30,
									 35, 116, 97, 112, 82),
						 draft1000=c(44, 25, 147, 127, 75,
									  95, 77, 118,95, 65,
									  57, 165, 140, 173, 119),
						  row.names=1897:1911)
