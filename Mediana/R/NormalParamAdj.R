######################################################################################################################

# Function: NormalParamAdj.
# Argument: p, Vector of p-values (1 x m)
#       par, List of procedure parameters: vector of hypothesis weights (1 x m) and correlation matrix (m x m)

# Description: Parametric multiple testing procedure based on a multivariate normal distribution

NormalParamAdj = function(p, par) {

  # Determine the function call, either to generate the p-value or to return description
  call = (par[[1]] == "Description")


  if (any(call == FALSE) | any(is.na(call))) {
    # Number of p-values
    p = unlist(p)
    m = length(p)


    # Extract the vector of hypothesis weights (1 x m) and correlation matrix (m x m)
    # If the first parameter is a matrix and no weights are provided, the hypotheses are assumed to be equally weighted
    if (is.null(par[[2]]$weight)) {
      w = rep(1/m, m)
    } else {
      w = unlist(par[[2]]$weight)
      if (is.null(par[[2]]$corr)) stop("Analysis model: Parametric multiple testing procedure: Correlation matrix must be specified.")
      corr = par[[2]]$corr
    }

    # Error checks
    if (length(w) != m) stop("Analysis model: Parametric multiple testing procedure: Length of the weight vector must be equal to the number of hypotheses.")
    if (sum(w) != 1) stop("Analysis model: Parametric multiple testing procedure: Hypothesis weights must add up to 1.")
    if (any(w < 0)) stop("Analysis model: Parametric multiple testing procedure: Hypothesis weights must be greater than 0.")

    if (sum(dim(corr) == c(m, m)) != 2) stop("Analysis model: Parametric multiple testing procedure: Correlation matrix is not correctly defined.")
    if (det(corr) <= 0) stop("Analysis model: Parametric multiple testing procedure: Correlation matrix must be positive definite.")


    # Compute test statistics based on a normal distribution
    stat = stats::qnorm(1 - p)

    # Adjusted p-values computed using a multivariate normal distribution function
    adjpvalue = sapply(stat, NormalParamDist, w, corr)


    result = adjpvalue
  }
  else if (call == TRUE) {
    if (is.null(par[[2]]$weight)) {
      w = rep(1/m, m)
    } else {
      w = unlist(par[[2]]$weight)
      if (is.null(par[[2]]$corr)) stop("Analysis model: Parametric multiple testing procedure: Correlation matrix must be specified.")
      corr = par[[2]]$corr
    }
    weight = paste0("Hypothesis weights={", paste(round(w, 3), collapse = ","),"}")
    corr = paste0("Correlation matrix={", paste(as.vector(t(corr)), collapse = ","),"}")
    result=list(list("Normal parametric multiple testing procedure"), list(weight,corr))
  }


  return(result)
}
# End of NormalParamAdj