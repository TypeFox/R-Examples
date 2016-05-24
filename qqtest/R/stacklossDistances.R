#' Mahalanobis squared distances of Brownlee's stack loss plant operation data based only on the explanatory variates (air flow, water temperature, and acid concentration).
#'
#' Mahalanobis distances were calculated using the mahalanobis R function. 
#' Under standard normal theory, these are approximately Chi-squared on 3 degrees of freedom.
#'
#' 
#' \code{with(stacklossDistances, qqtest(robust,dist="chi", df=3))} will show "outliers".
#' 
#' \code{with(stacklossDistances, qqtest(ordinary,dist="chi", df=3))} will show "inliers".
#'
#'
#' @format A data frame with 21 rows and 2 variates:
#' \describe{
#'   \item{ordinary}{Mahalanobis squared distances from the arithemetic mean using the eliptical contours of the sample covariance matrix.}
#'   \item{robust}{As with distances, these are Mahalobis squared distances butnow based on robust measures of location and covariance matrix (as determined from the default covRob of the robust package).}
#' }
#' @source
#' "Statistical Theory and Methodology in Science and Engineering", 
#' K.A. Brownlee, (1960, 2nd ed. 1965), Wiley, New York pp. 491-500.  
"stacklossDistances"

stacklossDistances <- data.frame(ordinary=c(5.07872843, 5.40443822, 2.53991916, 1.61772391, 0.09202564,
											0.59736629, 3.43235424, 3.43235424, 1.85129174, 3.04850429,
											2.14828274, 3.39113835, 2.19824823, 3.16407660, 2.85691631,
											1.66909308, 7.29008900, 2.25947354, 2.53835163, 0.65133605,
											4.73828830),
						   		robust=c(17.4057053, 18.4700184,  9.6049314,  0.9674260,  0.5768629,
						   						   0.6180002,  1.3463439,  1.3463439,  0.5884851,  2.7562423,
						   						   2.0204779,  2.8374068,  1.8862127,  3.3058572,  4.2819184,
						   						   2.4623482,  6.5266003,  2.0402159,  2.3659739,  0.5659061,
						   						   8.2772591),
						  		row.names=1:21
						    	)
