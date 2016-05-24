is.acct <- function(x){ # test for acc.fct
  if(!is.character(x)){stop("The acc.fct you entered is not of type character!")}
  if(!is.na(x[2])){stop("The acc.fct you entered is a vector!")}
  if(!(x %in% c("logistic") )){stop("The acc.fct you entered is not valid!")}
}