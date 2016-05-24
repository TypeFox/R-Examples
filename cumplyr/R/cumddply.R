cumddply <- function(data, equality.variables, inequality.variables, func)
{
  # equality.variables and inequality.variables must be disjoint.
  all.variables <- c(equality.variables, inequality.variables)
  
  # For all constraining variables, calculate the unique, sorted values of that variable.
  # Place these into an environment.
  local.env <- new.env()
  for (variable in all.variables)
  {
    assign(variable, sort(unique(get(variable, data))), envir = local.env)
  }
  
  # Find Cartesian product of all unqiue, sorted values of all variables
  cartesian.product <- cartesian_product(all.variables, envir = local.env)
  
  results <- data.frame()
  
  # Iterate over elements of the Cartesian product
  for (row.index in 1:nrow(cartesian.product))
  {
    local.data <- data
    
    for (variable in equality.variables)
    {
      local.data <- subset(local.data, get(variable, local.data) == cartesian.product[row.index, variable])
    }
    
    for (variable in inequality.variables)
    {
      local.data <- subset(local.data, get(variable, local.data) <= cartesian.product[row.index, variable])
    }
    
    results <- rbind(results, cbind(cartesian.product[row.index, ], func(local.data)))
    # Add row of new data to results containging value of function and current element from Cartesian product    
  }
  
  return(results)
}
