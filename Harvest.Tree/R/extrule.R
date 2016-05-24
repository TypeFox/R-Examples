#' Bound of rules
#' 
#' This function takes in a ruleset and output the lower and upper bounds of each rule.
#' @param myrules A 3 column matrix output of function "hughs.path.rpart"
#' @param varname the names of x variables
#' @return A p*2 matrix, p is the length of varname. The first column is the lower bound, the second column is the upper bound. The default lower bound is "-Inf",the default upper bound is "Inf". row corresponse to x variables ordered in the data matrix given to rpart.




extrule <- function(myrules, varname) 
{ 
  varnum <- length(varname)
  exrule <- data.frame(rep(-Inf,varnum),rep(Inf,varnum),rep(NA,varnum))
  for (i in 1:varnum) 
  { 
    for (j in nrow(myrules):1) 
    {
      if (myrules[j,1]==varname[i] & myrules[j,2]==0) {
        exrule[i,2] <- as.numeric(myrules[j,3])
        break }              
    }
    for (j in nrow(myrules):1) 
    {
      if (myrules[j,1]==varname[i] & myrules[j,2]==1) {
        exrule[i,1] <- as.numeric(myrules[j,3])
        break }              
    }  
    for (j in nrow(myrules):1) 
    {
      if (myrules[j,1]==varname[i] & myrules[j,2]==2) {
        exrule[i,3] <- myrules[j,3]
        break }              
    }            
  } 
  return(exrule)
}



