tree2glm <- function(xtree, xdata, Y.name, X.names, family = "binomial")
{
  if(!inherits(xtree, 'rpart')) stop('xtree have to be an rpart object!')
  indicators_tree <- sapply(tree2indicators(xtree), function(u) return(paste("as.integer(", u, ")")))
  
  nber_indicators <- length(indicators_tree)
  
  if(nber_indicators > 1) xformula <- as.formula(paste(Y.name, "~", paste(X.names, collapse = " + "), "+", paste(indicators_tree[-nber_indicators], collapse = "+")))
  else xformula <- as.formula(paste(Y.name, "~", paste(X.names, collapse = " + ")))
  
  fit <- glm(xformula, data = xdata, family = family)
  
  return(fit)
}
