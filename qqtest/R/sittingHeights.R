#' Sitting height in inches of female adults (aged 23-50).  
#'
#' From measurements made by Francis Galton at London's International Health Exhibition in 1884, published in 1885.
#'
#' 
#' 
#' \code{with(sittingHeights, plot(percentCumulative,upperBound, type="o", lwd=2, xlim=c(0,100), ylim=c(20,40),xlab="Percentage of women", ylab="Sitting heights in inches"))} will effect Galton's Ogive.
#' 
#' \code{with(sittingHeights, qqtest(binCentre, dist="normal", p = proportionCumulativeAdjusted, np=775, main="Sitting heights of women in inches"))} will effect a normal qqplot for this data.
#'
#'
#' @format A data frame with 9 rows and 8 variates, the first 5 of which are as recorded by Galton:
#' \describe{
#'   \item{lowerBound}{Sitting height greater than or equal to this lower bound in inches.}
#'   \item{upperBound}{Sitting height strictly less than this upper bound in inches.}
#'   \item{nCases}{Number of cases observed with sitting height between the two bounds.}
#'   \item{nCasesCumulative}{Number of cases observed with sitting height up to but not including the upper bound.}
#'   \item{percentCumulative}{Cumulative number of cases expressed as a percent.}
#'   \item{binCentre}{Average of the lower and upper bounds.}
#'   \item{nCasesCentred}{Number of cases assigned to the binCentre; half of cases observed with sitting height between the two bounds is assigned to be below the binCentre, half above (avoids producing 100 percent for last entry).}
#'   \item{proportionCumulativeAdjusted}{Cumulative proportions using nCasesCentred. }
#' }
#' @source
#' "The Application of a Graphic Method to Fallible Measures", 
#' Francis Galton, (1885), Journal of the Statistical Society of London, Jubilee Volume (June 22-24, 1885), pp. 262-265.  
"sittingHeights"

sittingHeights <- data.frame(lowerBound=c(29, 30, 31, 32, 33, 34, 35, 36, 37),
						   		upperBound=c(30, 31, 32, 33, 34, 35, 36, 37, 38),
						   		nCases=c(2, 8, 52, 116, 226, 227, 108, 31, 5),
						   		nCasesCumulative=c(2, 10, 62, 178, 404, 631, 739, 770, 775),
						   		percentCumulative=c(0.02, 1.3, 8.0, 23.0, 52.2, 81.4, 95.3, 99.4, 100.0),
								binCentre=c(29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5),
						   		nCasesCentred=c(1, 1+4, 4+26, 26+58, 58+113, 113+113.5, 113.5+54, 54+15.5, 15.5+2.5),
						  		proportionCumulativeAdjusted =c(1, 6, 36, 120, 291, 517.5, 685, 754.5, 772.5)/775,
						  		row.names=c(
						    	"29 to under 30 inches", 
						    	"30 to under 31 inches", 
						   	"31 to under 32 inches", 
						    	"32 to under 33 inches", 
						    	"33 to under 34 inches", 
						    	"34 to under 35 inches", 
						    	"35 to under 36 inches", 
						    	"36 to under 37 inches", 
						    	"37 to under 38 inches")
						    	)
