#' Strength of pull for 519 males aged 23-26.  
#'
#' From measurements made by Francis Galton at the International Health Exhibition in 1884.
#'
#' 
#' 
#' \code{qqtest(pullstrength$strength, p=pullstrength$percentCumulative/100, np=519, dist="uniform", main="Galton's ogive of pull strength for 519 males aged 23-26", xlab="Cumulative Proportions (Adjusted)",yAxisAsProbs=FALSE, ylab="Strength in lbs.", type="o")} will effect Galton's Ogive.
#' 
#' \code{qqtest(pullstrength$strength, p=pullstrength$percentAdjustedCumulative/100, np=519, dist="normal", main="Gaussian qqplot of pull strength for 519 males aged 23-26", xlab="Cumulative Proportions (Adjusted)",yAxisAsProbs=FALSE, ylab="Strength in lbs.", type="o")} will effect a normal qqplot for this data.
#'
#'
#' @format A data frame with 7 rows and 4 variates:
#' \describe{
#'   \item{strength}{Pull strength lower bound in pounds.}
#'   \item{nCases}{Number of cases observed with pull strength between this bound and the next.}
#'   \item{percentCases}{Percent of cases observed with pull strength between this bound and the next.}
#'   \item{percentCumulative}{Cumulative percent of cases observed with pull strength up to this bound.}
#'   \item{percentAdjustedCumulative}{Adjust Galton's cumulative percent to include only half  the cases between this bound and the next.}
#' }
#' @source
#' "Natural Inheritance", 
#' Francis Galton, (1889), Table 1, page 199.  
"pullstrength"

pullstrength <- data.frame(strength=c(40, 50, 60, 70, 80, 90, 100),
						   nCases=c(10, 42, 140, 168, 113,  22, 24),
						   percentCases=c(2, 8, 27, 33, 21, 4, 5),
						   percentCumulative=c(2, 10, 37, 70, 91, 95, 100),
						   percentAdjustedCumulative=c(1, 6, 23.5, 53.5, 80.5, 93, 97.5),
						   row.names=c(
						    "Under 50 lbs.", 
						    "Under 60 lbs.", 
						    "Under 70 lbs.", 
						    "Under 80 lbs.",
						    "Under 90 lbs.",
						    "Under 100 lbs.",
						    "Above 100 lbs.")
						    )
