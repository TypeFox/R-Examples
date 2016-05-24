#' CHDAGE data
#' 
#' chdage dataset.
#' 
#' @format A data.frame with 100 rows and 4 variables:
#' \describe{
#' \item{id}{Identification code (1 - 100)}
#' \item{age}{Age (Years)}
#' \item{agegrp}{Age group (1: 20-39, 2: 30-34, 3: 35-39, 4: 40-44, 5:
#' 45-49, 6: 50-54, 7: 55-59, 8: 60-69)}  
#' \item{chd}{Presence of CHD (1: No, 2: Yes)}
#' }
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(chdage,  n = 10)
#' summary(chdage)
#'
#' ## Figure 1.1 p. 5
#' plot(as.integer(chd)-1 ~ age,
#'      pch = 20,
#'      main = "Figure 1.1 p. 5",
#'      ylab = "Coronary heart disease",
#'      xlab = "Age (years)",
#'      data = chdage)
#'
#' ## Table 1.2
#' with(chdage, addmargins(table(agegrp)))
#' with(chdage, addmargins(table(agegrp, chd)))
#' (Means <- with(chdage, tapply(as.integer(chd)-1, list(agegrp), mean)))
#' 
#' ## Figure 1.2 p. 6
#' midPoints <- c(24.5, seq(32, 57, 5), 64.5)
#' plot(midPoints, Means, pch = 20,
#'      ylab = "Coronary heart disease (mean)",
#'      xlab = "Age (years)", ylim = 0:1,
#'      main = "Figure 1.2 p. 6")
#' lines(midPoints, Means)
#' 
#' ## Table 1.3
#' summary( mod1.3 <- glm( chd ~ age, family = binomial, data = chdage ))
#'
#' ## Table 1.4
#' vcov(mod1.3)
#'
#' ## Computing OddsRatio and confidence intervals for age ...
#' exp(coef(mod1.3))[-1]
#' exp(confint(mod1.3))[-1, ]
"chdage"
