
validate.covar <- function(null, resp.var){
  
  # give a warning is some levels of factor have small sample size
  check.small.level(null, resp.var)
  
  # remove dummy variables with SD = 0
  # Note that we need to excluded samples with missing at covariates
  # so some covariates having different values may become constant after removing some samples
  # we therefore update data frame 'null' by removing those columns
  null <- remove.const.covar(null, resp.var)
  
  # check multicollinearity
  null <- remove.redendant.covar(null, resp.var)
  
  null
  
}

