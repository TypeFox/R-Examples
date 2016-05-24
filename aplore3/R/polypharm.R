#' POLYPHARM data
#' 
#' polypharm dataset.
#' 
#' @format A data.frame with 3500 rows and 14 variables:
#' \describe{
#' \item{id}{Subject ID (1 - 500)}
#' \item{polypharmacy}{Outcome; taking drugs from more than three different
#' classes (1: No, 2: Yes)}
#' \item{mhv4}{Number of outpatient Mental Health Visits (1: none, 2: one
#' to five, 3: six to fourteen, 4: greater than 14)}
#' \item{inptmhv3}{Number of inpatient Mental Health Visits (1: none, 2:
#' one, 3: more than one)}
#' \item{year}{Year (2002 to 2008)}
#' \item{group}{Group (1: Covered Families and Children - CFC, 2: Aged,
#' Blind or Disabled - ABD, 3: Foster Care - FOS)}
#' \item{urban}{Location (1: Urban, 2: Rural)}
#' \item{comorbid}{Comorbidity (1: No, 2: Yes)}
#' \item{anyprim}{Any primary diagnosis (bipolar, depression, etc.) (1:
#' No, 2: Yes)}
#' \item{numprim}{Number of primary diagnosis (1: none, 2: one, 3: more
#' than one)}
#' \item{gender}{Gender (1: Female, 2: Male)}
#' \item{race}{Race (1: White, 2: Black, 3: Other)}
#' \item{ethnic}{Ethnic category (1: Non-Hispanic, 2: Hispanic)}
#' \item{age}{Age (Years and months, two decimal places)}
#' }
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(polypharm, n = 10)
#' summary(polypharm)
"polypharm"
