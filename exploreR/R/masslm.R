#' Mass Linear Regression
#'
#' This function takes in a dataframe, the dependent variable, and optionally a character vector of independent variables you want the function to ignore.  It then produces a dataframe of regression results. 
#'
#' @param data data.frame object that contains both the dependent variable and predictor variables you want to regress.
#' 
#' @param dv.var single dependent variable you want to regress your predictors on. 
#' 
#' @param ignore accepts a character vector of one or more variables you want the function to skip.  If nothing is passed through this option, the function will attempt to run a regression between the dependent variable and every other column of data.  
#' 
#' @param p.round set to TRUE by default.  If left TRUE, will round off the P.value outputs to their 6 significant digits.  Can be a problem for numbers larger than 999999, set to false to return the raw number.
#'
#' @param c.round set to TRUE by default.  If left TRUE, will round off the Coefficient outputs to their 6 significant digits.  Can be a problem for numbers larger than 999999, set to false to return the raw number.
#'
#' @return data.frame containing three columns of data, IV, Coefficient, and P.Value.  If one of the columns of data not excluded from the function contained character type data, the function will print an error recommending the user attempt to convert the variable to a factor.
#'
#' @export
#' 
#' @examples 
#' exam.df <- iris
#' masslm(exam.df, "Sepal.Width", ignore = "Species")
#' masslm(exam.df, "Sepal.Width", ignore = c("Species", "Petal.Width"))

masslm <- function (data, dv.var, ignore = NULL, p.round = TRUE, c.round = TRUE) {
  
  # Things that happen before the function takes place
  mydata <- data
  DV <- dv.var
  vars.list <- colnames(mydata)
  models <- data.frame(IV=character(),
                       Coefficient=double(),
                       P.value=double(),
                       R.squared=double())  
  
  # Searches through the list of vars and compares them to the list of ignores and the dv.  The function skips those.
  for (i in vars.list) {
    if(is.character(mydata[,i])) { 
      warning(i, " appears to be character type.  Try converting to factor using as.factor(mydata$", i, ")")
    } else {
      if (!(i %in% ignore | i == dv.var )) {
        IV <- i
        Coefficient <- summary(stats::lm(paste(dv.var, " ~ ", i, sep = ""), data = mydata))$coefficients[2, "Estimate"]
        P.value <- summary(stats::lm(paste(dv.var, " ~ ", i, sep = ""), data = mydata))$coefficients[2, "Pr(>|t|)"]
        R.squared <- summary(stats::lm(paste(dv.var, " ~ ", i, sep = ""), data = mydata))$r.squared
        models <- rbind(models, cbind(IV, Coefficient, P.value, R.squared))
      }
    }
  } 
  models$IV <- as.character(models$IV)
  models$P.value <- as.numeric(as.character(models$P.value))
  models$Coefficient <- as.numeric(as.character(models$Coefficient))
  models$R.squared <- as.numeric(as.character(models$R.squared))
  if(p.round) models$P.value <- signif(models$P.value, 4)
  if(c.round) models$Coefficient <- signif(models$Coefficient, 4)
  return(models)
}