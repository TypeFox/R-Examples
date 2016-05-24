#' GLOW500 data
#' 
#' glow500 dataset.
#' 
#' @format A data.frame with 500 rows and 15 variables:
#' \describe{
#' \item{sub_id}{Identification Code (1 - n)}
#' \item{site_id}{Study Site (1 - 6)}
#' \item{phy_id}{Physician ID code (128 unique codes)}
#' \item{priorfrac}{History of Prior Fracture (1: No, 2: Yes)}
#' \item{age}{Age at Enrollment (Years)}
#' \item{weight}{Weight at enrollment (Kilograms)}
#' \item{height}{Height at enrollment (Centimeters)}
#' \item{bmi}{Body Mass Index (Kg/m^2)}
#' \item{premeno}{Menopause before age 45 (1: No, 2: Yes)}
#' \item{momfrac}{Mother had hip fracture (1: No, 2: Yes)}
#' \item{armassist}{Arms are needed to stand from a chair (1: No, 2: Yes)}
#' \item{smoke}{Former or current smoker (1: No, 2: Yes)}
#' \item{raterisk}{Self-reported risk of fracture (1: Less than others of
#' the same age, 2: Same as others of the same age, 3: Greater than others
#' of the same age)}
#' \item{fracscore}{Fracture Risk Score (Composite Risk Score)}
#' \item{fracture}{Any fracture in first year (1: No, 2: Yes)}
#' }
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(glow500, n = 10)
#' summary(glow500)
#' 
#' ## Table 2.2 p. 39
#' summary(mod2.2 <- glm(fracture ~ age + weight + priorfrac +
#'                                  premeno + raterisk,
#'                       family = binomial,
#'                       data = glow500))
#'
#' ## Table 2.3 p. 40
#' summary(mod2.3 <- update(mod2.2, . ~ . - weight - premeno))
#'
#' ## Table 2.4 p. 44
#' vcov(mod2.3)
#'
#' ## Table 3.6 p. 58
#' contrasts(glow500$raterisk)
#'
#' ## Contrasts: Table 3.8 and 3.9 p. 60
#' contrasts(glow500$raterisk) <- matrix(c(-1,-1,1,0,0,1), byrow= TRUE, ncol = 2)
#' summary(mod3.9 <- glm(fracture ~ raterisk, family = binomial,
#'                       data = glow500))
#' # cleaning modified dataset ...
#' rm(glow500)
#'
#' ## Table 5.1 pg 160 - Hosmer-Lemeshow test (with vcdExtra package)
#' mod4.16 <- glm(fracture ~ age * priorfrac + height + momfrac * armassist +
#'                           I(as.integer(raterisk) == 3) ,
#'                family = binomial,
#'                data = glow500)
#' library(vcdExtra)
#' summary(HLtest(mod4.16))
#'
#' ## Table 5.3 p. 171 - Classification table
#' glow500$pred4.16 <- predict(mod4.16, type = "response")
#' with(glow500, addmargins(table( pred4.16 > 0.5, fracture)))
#'
#' ## Sensitivy, specificity, ROC (using pROC)
#' library(pROC)
#'
#' ## Figure 5.3 p. 177 - ROC curve (using pROC package)
#' print(roc4.16 <- roc(fracture ~ pred4.16, data = glow500))
#' plot(roc4.16, main = "Figure 5.3 p. 177")
#'
#' ## Table 5.8 p. 175
#' vars <- c("thresholds","sensitivities","specificities")
#' tab5.8 <- data.frame(roc4.16[vars])
#' ## Now, for printing/comparison purposes, steps below in order to find
#' ## threshold values most similar to those in the table
#' findIndex <- function(x, y) which.min( (x-y)^2 )
#' cutPoints <- seq(0.05, 0.75, by = 0.05)
#' tableIndex <- mapply(findIndex, y = cutPoints,
#'                      MoreArgs = list(x = roc4.16$thresholds))
#' ## And finally, let's print a reasonable approximation of table 5.8
#' writeLines("\nTable 5.8 p. 175\n")
#' tab5.8[tableIndex, ]
#'
#' ## Figure 5.1 p. 175
#' plot(specificities ~ thresholds, xlim = c(0, 1), type = "l",
#'      xlab = "Probabilty cutoff", ylab = "Sensitivity/specificity",
#'      ylim = c(0, 1), data = tab5.8, main = "Figure 5.1 p. 175")
#' with(tab5.8, lines(thresholds, sensitivities, col = "red"))
#' legend(x = 0.75, y = 0.55, legend = c("Sensitivity", "Specificity"),
#'        lty = 1, col = c("red","black"))
#' abline(h = c(0, 1), col = "grey80", lty = "dotted")
"glow500"

#' GLOW11M data
#' 
#' glow11m dataset.
#' 
#' @format A data.frame with 238 rows and 16 variables: the covariate are
#' the same as those from \code{\link{glow500}} with the addition of
#' \describe{
#' \item{pair}{Pair Identification Code (1-119)}
#' }
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(glow11m, n = 10)
#' summary(glow11m)
#'
#' ## Table 7.2 p. 252
#' library(survival)
#' mod7.2 <- clogit(as.numeric(fracture) ~ height + weight + bmi +
#'                  priorfrac + premeno + momfrac + armassist + raterisk +
#'                  strata(pair), data = glow11m)
#' summary(mod7.2)
"glow11m"


#' GLOW_BONEMED data
#' 
#' glow_bonemed dataset.
#' 
#' @format A data.frame with 500 rows and 18 variables: the covariate are
#' the same as those from \code{\link{glow500}} with the addition of
#' \describe{
#' \item{bonemed}{Bone medications at enrollment (1: No, 2: Yes)}
#' \item{bonemed_fu}{Bone medications at follow-up (1: No, 2: Yes)}
#' \item{bonetreat}{Bone medications both at enrollment and follow-up (1: No,
#' 2: Yes)}
#' }
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(glow_bonemed, n = 10)
#' summary(glow_bonemed)
"glow_bonemed"


#' GLOW_MIS_COMP data
#' 
#' glow_mis_comp dataset.
#' 
#' @format A data.frame with 500 rows and 10 variables: the covariate are
#' the same as those from \code{\link{glow500}}, without \code{bmi},
#' \code{premeno}, \code{armassist}, \code{smoke} and \code{fracscore}.
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(glow_mis_comp, n = 10)
#' summary(glow_mis_comp)
"glow_mis_comp"

#' GLOW_MIS_WMISSING data
#' 
#' glow_mis_wmissing dataset.
#' 
#' @format A data.frame with 500 rows and 10 variables: the covariate are
#' the same as those from \code{\link{glow500}}, without \code{bmi},
#' \code{premeno}, \code{armassist}, \code{smoke} and \code{fracscore}.
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(glow_mis_wmissing, n = 10)
#' summary(glow_mis_wmissing)
"glow_mis_wmissing"

#' GLOW_RAND data
#' 
#' glow_rand dataset.
#' 
#' @format A data.frame with 500 rows and 15 variables: the covariate are
#' the same as those from \code{\link{glow500}}.
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#' Logistic Regression, 3rd ed., New York: Wiley
#' @examples
#' head(glow_rand, n = 10)
#' summary(glow_rand)
"glow_rand"
