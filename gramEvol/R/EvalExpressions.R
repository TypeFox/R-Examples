EvalExpressions <- function(expressions, envir = parent.frame()) {
  # evaluates one or a collection of expression to numerial value
  
  vals <- list()

  # handle the expression collections
  if (class(expressions) == "expression") {
  
    expressions = as.list(expressions)
    for (expr in expressions){
      result <- eval(expr, envir=envir)
      vals <- c(vals, list(result))
    }
  } else {
    # if the input is an expression, make it in to a list
    if (class(expressions) == "GEPhenotype") {
      expressions = as.list(expressions)
    }

    # evaluate the expressions
    vals <- list()
    for (expr in expressions){
      
      # either use phenotype class
      if (class(expr) == "GEPhenotype") {
        if (expr$type == "T") {
          # extract the expression for terminal grammar
          expr = expr$parsed
        } else {
          warning("Can not evaluate non-terminal expression")
          expr = "NA"
        }
        # or evaluate character string
      } else if (class(expr) == "character") {
        expr = parse(text = expr)
      } else if (class(expr) != "expression") {
        warning("Can not evaluate invalid expression")
        expr = "NA"
      }
      
      result <- eval(expr, envir=envir)
      vals <- c(vals, list(result))
    }
  }

  # return a vector if a single expression
  if (length(expressions) == 1)
    return (vals[[1]])
  
  # convert list to dataframe and set column names
  vals = do.call(data.frame, vals)
  colnames(vals) = paste0("expr", seq_along(expressions))
  
  return (vals)
}
