smc <-
function(data)
{
  are.na = FALSE
  coeffs = data.frame(matrix(ncol = ncol(data)))
  
  for(i_dep in names(data))
  {
    dummy.data.frame = data.frame(matrix(ncol = length(names(data)), nrow = 1))
    colnames(dummy.data.frame) = c(i_dep, names(data)[-which(names(data) == i_dep)])
    fm = formula(dummy.data.frame)
    fit = lm(fm, data = data)
    if(any(is.na(fit$coefficients))){
      corr.var = names(fit$coefficients[is.na(fit$coefficients)])
      stop(paste(corr.var, "is perfectly correlated with a prevois measurment. It must be removed."))
    }
    coeffs[i_dep,] = fit$coefficients
  }
  coeffs = coeffs[-1,]
  if(are.na)
  {
      print("**ATTENTION**")
      print("You have perfectly correlated data among your measurments.")
      print("We recommend that you discard one of them.")
      print("Use 'check.correlation()' to find the correlated variables")
  }
  
  mult.corr.coeff = NULL
  
  for(i_dep in rownames(coeffs))
  {
    r = mcc(data[,-which(names(data) == i_dep)], data[,which(names(data) == i_dep)], coeffs[i_dep,]  )
    
    mult.corr.coeff = c(mult.corr.coeff, r)
  }
  
  return(mult.corr.coeff)
  
}
