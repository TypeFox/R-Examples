#' @import TeachingSampling
#' @export
#' 
#' @title
#' Estimated sample Effects of Design (DEFF)
#' @description 
#' This function returns the estimated design effects for a set of inclusion probabilities and the variables of interest. 
#' @details
#' The design effect is somehow defined to be the ratio between the variance of a complex design and the variance of a simple design. 
#' When the design is stratified and the allocation is proportional, this measures reduces to
#' \deqn{DEFF_{Kish} = 1 + CV(w)}
#' where w is the set of weights (defined as the inverse of the inclusion probabilities) along the sample, and CV refers 
#' to the classical coefficient of variation. Although this measure is #' motivated by a stratified sampling design, 
#' it is commonly applied to any kind of survey where sampling weight are unequal. On the other hand, the Spencer's DEFF 
#' is motivated by the idea that a set of weights may be efficent even when they vary, and is defined by:
#' \deqn{DEFF_{Spencer} = (1 - R^2) * DEFF_{Kish} + \frac{\hat{a}^2}{\hat{\sigma}^2_y} * (DEFF_{Kish} - 1) }
#' where 
#' \deqn{\hat{\sigma}^2_y = \frac{\sum_s w_k (y_k - \bar{y}_w)^2}{\sum_s w_k}}
#' and \eqn{\hat{a}} is the estimation of the intercept in the following model
#' \deqn{y_k = a + b * p_k + e_k}
#' with \eqn{p_k = \pi_k / n} is an standardized sampling weight. Finnaly, \eqn{R^2} is the R-squared of this model. 
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param y Vector, matrix or data frame containing the recollected information of the variables of interest for every unit in the selected sample.
#' @param pik Vector of inclusion probabilities for each unit in the selected sample.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas.
#' Valliant, R, et. al. (2013), \emph{Practical tools for Design and Weighting Survey Samples}. Springer
#' @examples 
#' #############################
#' # Example with BigLucy data #
#' #############################
#' data(BigLucy)
#' attach(BigLucy)
#' 
#' # The sample size
#' n <- 400
#' res <- S.piPS(n, Income)
#' sam <- res[,1]
#' # The information about the units in the sample is stored in an object called data
#' data <- BigLucy[sam,]
#' attach(data)
#' names(data)
#' # Pik.s is the inclusion probability of every single unit in the selected sample
#' pik <- res[,2]
#' # The variables of interest are: Income, Employees and Taxes
#' # This information is stored in a data frame called estima
#' estima <- data.frame(Income, Employees, Taxes)
#' E.piPS(estima,pik)
#' DEFF(estima,pik)

DEFF <- function(y, pik){
  y <- cbind(1,y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  DEFFs <- matrix(NA,nrow=dim(y)[2],ncol=2)
  colnames(DEFFs)=c("DEFF.Kish", "DEFF.Spencer")
  rownames(DEFFs)<-names(y)
  
  n <- length(pik)
  wk = 1/pik
  pk = pik/n
  DEFF.Kish = 1 + var(wk) / mean(wk)^2
  
  for(k in 1:dim(y)[2]){
  model.pk <- lm(y[,k] ~ pk)
  R2 <- summary(model.pk)$r.squared
  a <- model.pk$coefficients[1]
  vary <- sum(wk * (y[,k] - (sum(y[,k])/sum(wk)))^2) / sum(wk)
  DEFF.Spen <- (1 - R2) * DEFF.Kish + (a^2 / vary) * (DEFF.Kish - 1)
  DEFFs[k,] <- c(DEFF.Kish, DEFF.Spen)
  }
  return(DEFFs) 
}