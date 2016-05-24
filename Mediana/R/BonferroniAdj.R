######################################################################################################################

# Function: BonferroniAdj.
# Argument: p, Vector of p-values (1 x m)
#           par, List of procedure parameters: vector of hypothesis weights (1 x m)
# Description: Bonferroni multiple testing procedure.

BonferroniAdj = function(p, par) {

  # Determine the function call, either to generate the p-value or to return description
  call = (par[[1]] == "Description")


  # Number of p-values
  m = length(p)

  # Extract the vector of hypothesis weights (1 x m)
  if (!any(is.na(par[[2]]))) {
    if (is.null(par[[2]]$weight)) stop("Analysis model: Bonferroni procedure: Hypothesis weights must be specified.")
    w = par[[2]]$weight
  } else {
    w = rep(1/m, m)
  }

  # Error checks
  if (length(w) != m) stop("Analysis model: Bonferroni procedure: Length of the weight vector must be equal to the number of hypotheses.")
  if (sum(w)!=1) stop("Analysis model: Bonferroni procedure: Hypothesis weights must add up to 1.")
  if (any(w < 0)) stop("Analysis model: Bonferroni procedure: Hypothesis weights must be greater than 0.")

  if (any(call == FALSE) | any(is.na(call))) {
    # Adjusted p-values
    adjpvalue = pmin(1, p/w)
    result = adjpvalue
  }
  else if (call == TRUE) {
    weight = paste0("Weight={",paste(round(w,2), collapse = ","),"}")
    result=list(list("Bonferroni procedure"),list(weight))
  }


  return(result)
}
# End of BonferroniAdj