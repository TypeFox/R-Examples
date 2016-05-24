

mgmfit <- function(
  data, # data matrix, col=variables
  type, # data type for col 1:ncol; c=categorical, g=gaussian, p=poisson, e=exponential
  lev, # number of categories of categorical variables, continuous variables have level=1
  lambda.sel = "EBIC", # method for penalization parameter (lambda) -selection 
  folds = 10, # folds in case CV is used for lambda selection
  gam = .25, # tuning parameter for EBIC, in case EBIC is used for lambda selection
  d = 2, # maximal degree of the true graph
  rule.reg = "AND", # parameter-aggregation of categorical variables
  pbar = TRUE, # shows a progress bar if TRUE
  method = 'glm',  # which method should be used for each nodewise regression?
  missings = 'error', # handling of missing data
  weights = NA, # weights for observations 
  ret.warn = TRUE # TRUE returns warnings, makes sense to switch off for time varying wrapper
)

{
  
  outlist <- mgmfit_core(data = data, 
                         type = type, 
                         lev = lev, 
                         lambda.sel = lambda.sel, 
                         folds = folds, 
                         gam = gam, 
                         d = d, 
                         rule.reg = rule.reg, 
                         pbar = pbar, 
                         method = method, 
                         missings = missings, 
                         weights = weights, 
                         ret.warn = ret.warn, 
                         VAR = FALSE) # use standard mgm.fit; no AR model
  
  return(outlist)
} 

