iddply <- function(data,
                   equality.variables = c(),
                   lower.bound.variables = c(),
                   upper.bound.variables = c(),
                   norm.ball.variables = list(),
                   func = function (df) {df})
{
  # All constraints should be disjoint.
  all.variables <- c(equality.variables,
                     lower.bound.variables,
                     upper.bound.variables,
                     names(norm.ball.variables))
    
  # For all constraining variables, calculate the unique, sorted values of that variable.
  # Place these into an environment.
  local.env <- new.env()
  for (variable in all.variables)
  {
    local.env[[variable]] <- sort(unique(get(variable, data)))
  }
  
  # Find Cartesian product of all unique, sorted values of all variables
  cartesian.product <- cartesian_product(all.variables, envir = local.env)
  
  # Iterate over elements of the Cartesian product
  results <- cartesian.product
  output <- func(data)
  output.col.size <- ncol(output)
  results <- cbind(results, output)
  
  for (row.index in 1:nrow(cartesian.product))
  {
    indices <- list()
    
    for (variable in equality.variables)
    {
      indices[[variable]] <- which(data[, variable] == cartesian.product[row.index, variable])
    }
    
    for (variable in lower.bound.variables)
    {
      indices[[variable]] <- which(data[, variable] >= cartesian.product[row.index, variable])
    }
    
    for (variable in upper.bound.variables)
    {
      indices[[variable]] <- which(data[, variable] <= cartesian.product[row.index, variable])
    }
    
    for (variable in names(norm.ball.variables))
    {
      indices[[variable]] <- which((data[, variable] - cartesian.product[row.index, variable])^2 <= norm.ball.variables[[variable]]^2)
    }
    
    index.set <- Reduce(intersect, indices)
    
    results[row.index, (ncol(cartesian.product) + 1):ncol(results)] <- func(data[index.set, ])
  }
  
  names(results) <- c(all.variables, paste('Var', seq_len(ncol(results) - length(all.variables)), sep = ''))
  return(results)
}
