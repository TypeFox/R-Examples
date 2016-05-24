QT <- function(model, significance = 0.05, hc=4, h0=0){
  if (class(model) != "lm") 
    stop("The argument model must have class lm.")
  if (significance >= 1 || significance <= 0) {
    stop("The significance level should belong to the open interval (0,1).")
  }
  if(length(grep("- 1",as.character(as.formula(model))))==1){
    stop("The model must contain the intercept.")
  }
  variables = labels(model)
  p = length(variables)
  n = length(model$residuals)
  statistics =  abs((as.numeric(model$coefficients) - h0)/sqrt(diag(HC(model,method=hc))))
  names(statistics) =  c("Intercept",variables)
  p_value = pnorm(statistics, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  result = list("model" = model, "statistics" = statistics, "p_value" = p_value)
  class(result) <- "list"
  return(result)
}