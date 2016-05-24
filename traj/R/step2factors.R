step2factors <-
function(trajMeasures, num.factors = NULL, discard  = NULL, verbose = TRUE, ...)
{

  # Deal with varibles to discard
  if(!is.null(discard)){
    if(class(discard) == "character")
      vars.to.discard = which(names(trajMeasures$measurments) %in% discard)
    else
      vars.to.discard = discard
    
    if(19 %in% vars.to.discard)
      stop("m18 will automatically be removed. Do not include it in the 'discard' variable.")
    if(length(vars.to.discard) != length(discard))
      
      stop("Not all variables in 'discard' are to be removed. There is an error in the format of 'discard'.")
    
    data = trajMeasures$measurments[,-vars.to.discard]
  }
  else
    data = trajMeasures$measurments
 
    
  
  # Sizing data
  dim.of.data = dim(data)
  sample.size = dim.of.data[1]
  
  # Deal with IDs
    IDvector = data[1] 
    data = data[-1]



  # Remove m18 id correslation larger than 0.95
  if(cor(data$m17,data$m18) >= 0.95){
    data = data[,-which(names(data) == "m18")]
  }
  
  # Deal with remaining correlated variables
  corr.vars = check.correlation(data, verbose = FALSE, is.return = TRUE)
  if(!is.null(corr.vars)){
    corr.vars.pos = which(names(data) %in% corr.vars[,1])
    data = data[,-corr.vars.pos]
    print(paste(corr.vars[,1], "is removed because it is perfectly correlated with", corr.vars[,2]))
  }

  # Checking validity of num.factors
  if(num.factors > ncol(data) && !is.null(num.factors))
    stop("Requesting more factors in 'num.factors' than available variables.")
  
  max.num.obs = dim(data)[2]
  
  eigen.values = NULL;
  pricipal.factors = NULL;
  
  # Calculate the number of factors to use
  if(is.null(num.factors))
  {
    if(verbose) 
      print("Computing reduced correlation e-values...")
    
    eigen.values = reduced.eigen(data)
    
    num.factors = length(which(eigen.values$values >= 1))
  }

  # Choose the principal varaibles that will represent the factors
  principal.factors = principal(data, rotate = "varimax", nfactors = num.factors, ...)
  principal.variables = c(rep(NA , num.factors))
  
  for(i_factors in 1: num.factors){
    principal.variables[i_factors] = which.max(abs(principal.factors$loadings[,i_factors]))
  }
  
  principal.variables = sort(principal.variables)
  
  # Bind the vectors of the factor variables to the ID vector
  output = IDvector

  for(i_col in 1 : num.factors){
    output = cbind(output, data[principal.variables[i_col]])
  }

  # Create list to export
  trajFactors = structure(list( factors  = output, e.values = eigen.values, princ.fact = principal.factors,
                  measurments = trajMeasures$measurments, data = trajMeasures$data, time = trajMeasures$time), class = "trajFactors")
  
  return(trajFactors)
}
