iddply2 <- function(data,
                   equality.variables = c(),
                   lower.bound.variables = c(),
                   upper.bound.variables = c(),
                   norm.ball.variables = list(),
                   func = function (df) {df})
{
  library('data.table')
  data <- data.table(data)

  # This approach is still not right, because I'm doing vector scans.
  # How do I get inequalities using data.table?
  
  # All constraints should be disjoint.
  all.variables <- c(equality.variables,
                     lower.bound.variables,
                     upper.bound.variables,
                     names(norm.ball.variables))
  
  setkeyv(data, all.variables)
    
  # For all constraining variables, calculate the unique, sorted values of that variable.
  # Place these into an environment.
  local.env <- new.env()
  for (variable in all.variables)
  {
    local.env[[variable]] <- sort(unique(get(variable, data)))
  }
  
  # Find Cartesian product of all unqiue, sorted values of all variables
  cartesian.product <- cartesian_product(all.variables, envir = local.env)
  
  # Iterate over elements of the Cartesian product
  results <- data.frame(stringsAsFactors = FALSE)
  
  for (row.index in 1:nrow(cartesian.product))
  {
    local.data <- data
    
    for (variable in equality.variables)
    {
      #local.data <- local.data[local.data[, variable] == cartesian.product[row.index, variable], ]
      # Generate expression from item above, substituting variable with its value
      # Need to use proper string literals for cartesian.product values.
      expr1 <- parse(text = paste(variable, '==', cartesian.product[row.index, variable]))
      local.data <- local.data[eval(expr1), ]
    }
    
    for (variable in lower.bound.variables)
    {
      #local.data <- local.data[local.data[, variable] >= cartesian.product[row.index, variable], ]    
      expr1 <- parse(text = paste(variable, '>=', cartesian.product[row.index, variable]))
      local.data <- local.data[eval(expr1), ]
    }
    
    for (variable in upper.bound.variables)
    {
      #local.data <- local.data[local.data[, variable] <= cartesian.product[row.index, variable], ]
      expr1 <- parse(text = paste(variable, '<=', cartesian.product[row.index, variable]))
      local.data <- local.data[eval(expr1), ]      
    }
    
    for (variable in names(norm.ball.variables))
    {
      local.data <- local.data[(local.data[, variable] - cartesian.product[row.index, variable])^2 <= norm.ball.variables[[variable]]^2, ]
      # Handle this one later.
    }
    
    results <- rbind(results, cbind(cartesian.product[row.index, ], func(local.data)))
  }
  
  names(results) <- c(all.variables, paste('Var', seq_len(ncol(results) - length(all.variables)), sep = ''))
  return(results)
}
