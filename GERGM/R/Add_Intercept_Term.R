add_intercept_term <- function(
  formula){
  formula<- as.formula(formula)
  parsed <- deparse(formula)
  if(length(parsed) > 1){
    parsed <- paste0(parsed, collapse = " ")
  }
  if(grepl("intercept",parsed)){
    stop("The intercept term is specified automatically, please remove the 'intercept' term from your formula.")
  }
  parsed <- paste(parsed," + intercept",sep = "")
  formula <- as.formula(parsed)
  return(formula)
}
