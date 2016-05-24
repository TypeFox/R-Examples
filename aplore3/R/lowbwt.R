#' LOWBWT data
#' 
#' lowbwt dataset.
#' 
#' @format A data.frame with 189 rows and 11 variables:
#' \describe{
#' \item{id}{Identification Code}
#' \item{low}{Low birth weight (1: >= 2500, 2: < 2500 g)}
#' \item{age}{Age of mother (Years)}
#' \item{lwt}{Weight of mother at last menstrual period (Pounds)}
#' \item{race}{Race (1: White, 2: Black, 3: Other)}
#' \item{smoke}{Smoking status during pregnancy (1: No, 2: Yes)}
#' \item{ptl}{History of premature labor (1: None, 2: One, 3: Two, etc)}
#' \item{ht}{History of hypertension (1: No, 2: Yes)}
#' \item{ui}{Presence of Uterine irritability (1: No, 2: Yes)}
#' \item{ftv}{Number of physician visits during the first trimester (1:
#' None, 2: One, 3: Two, etc)} 
#' \item{bwt}{Recorded birth weight (Grams)}
#' }
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(lowbwt, n = 10)
#' summary(lowbwt)
"lowbwt"
