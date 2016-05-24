#' Bacteria from Delaware River water entering the Torresdale Filter 
#' of the Philadelphia water supply 1913.  
#'
#' The number of bacteria per cubic centimetre was measured daily 
#' for river water entering the Torresdale filter of the Philadelphia
#' water supply through the whole of 1913.  
#' 
#' Rather than the individual daily results, recorded here are 22 values 
#' together with the percentage of days whose value was less than or equal
#' to the recorded value.
#' These quantiles are therefore based on 365 daily measurements.
#' 
#' Values were taken from a logarithmic-normal probability plot dated 
#' January 21, 1915 as it appeared in Figure 22 of George C. Whipple's 
#' 1916 (Part 2) paper on the ``Element of Chance in Sanitation''.
#' This paper introduces logarithimic-normal probability paper (or a log-normal qqplot).
#' 
#' \code{with(bacteria, plot(qnorm(percentTime/100), log(count,10), type="o"))} 
#' reproduces Whipple's 1915 plot.
#' 
#' \code{with(bacteria, qqtest(data=count, p=percentTime/100, np=365, dist="log-normal", type="o"))}
#'  will effect a qqtest plot for this data.  More detail can be had from:
#' \code{with(bacteria, qqtest(data=log(count, 10), p=percentTime/100, np=365, dist="normal", type="o"))}
#'
#'
#' @format A data frame with 22 rows and 2 variates:
#' \describe{
#'   \item{count}{Number of bacteria per cc.}
#'   \item{percentTime}{Percent of days (out of 365) having bacteria count 
#' less than or equal to the measured count.}
#' }
#' @source 
#' "The Element of Chance in Sanitation", 
#' George C. Whipple, Journal of the Franklin Institue, Volume 182,
#' July and August (1916), pp. 37-59 and 205-227.
#' Data taken directly from Figure 22, page 209.  
"bacteria"

bacteria <- data.frame( count = c(680, 720, 775, 810,875, 
			  						  2000, 2900, 3000, 5000,8000,
			  						  10000, 12000, 15000, 
			  						  20000, 30000,
			  						  40000, 41000, 41000, 
			  						  42000, 47000, 47000, 80500),
			  			   percentTime = c(0.12, 0.4, 0.7,0.75, 1, 
											  12.5, 26, 27, 51, 73, 
											  83.5, 84.5, 88.5,
											  91.5, 95.5, 
											  98.3, 98.5, 98.75,
											  99, 99.25, 99.6, 99.85 )
							)
