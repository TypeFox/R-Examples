## #########################################
## Function to estimate the propensity score
## using a logistic model
## #########################################

pscore <- function(formula,
                   data,
                   family=binomial,
                   na.action=na.exclude,
                   name.pscore="pscore",
                   ...
                   )
{
  ## Check data argument
  if (missing(data)){
    stop("Argument 'data' is missing.")
  }else{
    if(!inherits(data,"data.frame"))
      stop("Argument 'data' is not of class 'data.frame'.")
  }

  ## Check formula argument
  if (missing(formula)){
    stop("Argument 'formula' is missing.")
  }else{
    if(!inherits(formula,"formula"))
      stop("Argument 'formula' is not of class 'formula'.")
  }

  ## Check arguments if necessary
  if(any(names(data)==name.pscore))
    stop(paste("Variable 'name.pscore'=",name.pscore," already exists in data.", sep=""))
  
  ## Extract treatment ( == response here)
  name.treat <- names(model.frame(formula,data))[1] 
  treat      <- data[,name.treat]

  ## Fit PS model
  ps.model <- glm(formula, family, data, na.action=na.action, ...)

  ## Predict PS values
  data[,name.pscore] <- predict(ps.model, type="response")

  ## Define output:
  pscore <- data[,name.pscore]
 
  output <- list(data           = data,
                 pscore         = pscore,
                 name.pscore    = name.pscore,
                 formula.pscore = formula,
                 model.pscore   = ps.model,
                 treat          = treat,
                 name.treat     = name.treat)
  
  class(output) <- c("pscore") 

  return(output)

}
