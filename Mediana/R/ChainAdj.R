######################################################################################################################

# Function: ChainAdj.
# Argument: p, Vector of p-values (1 x m)
#           par, List of procedure parameters: vector of hypothesis weights (1 x m) matrix of transition parameters (m x m)
# Description: Chain multiple testing procedure.

ChainAdj = function(p, par) {

  # Determine the function call, either to generate the p-value or to return description
  call = (par[[1]] == "Description")

  if (any(call == FALSE) | any(is.na(call))) {
    # Number of p-values
    m = length(p)
    # Extract the vector of hypothesis weights (1 x m) and matrix of transition parameters (m x m)
    if (is.null(par[[2]]$weight)) stop("Analysis model: Chain procedure: Hypothesis weights must be specified.")
    if (is.null(par[[2]]$transition)) stop("Analysis model: Chain procedure: Transition matrix must be specified.")
    w = par[[2]]$weight
    g = par[[2]]$transition

    # Error checks
    if (sum(w)!=1) stop("Analysis model: Chain procedure: Hypothesis weights must add up to 1.")
    if (any(w<0)) stop("Analysis model: Chain procedure: Weights must be greater than 1.")
    if (sum(dim(g) == c(m, m)) != 2)
      stop("Analysis model: Chain procedure: The dimension of the transition matrix is not correct.")
    if (any(rowSums(g)>1))
      stop("Analysis model: Chain procedure: The sum of each row of the transition matrix must be lower than 1.")
    if (any(g < 0))
      stop("Analysis model: Chain procedure: The transition matrix must include only positive values.")

    pmax = 0


    # Index set of processed, eg, rejected, null hypotheses (no processed hypotheses at the beginning of the algorithm)
    processed = rep(0, m)


    # Adjusted p-values
    adjpvalue = rep(0, m)


    # Loop over all null hypotheses
    for (i in 1:m) {
      # Find the index of the smallest weighted p-value among the non-processed null hypotheses
      ind = argmin(p, w, processed)
      if (ind>0){
        adjpvalue[ind] = max(p[ind]/w[ind], pmax)
        adjpvalue[ind] = min(1, adjpvalue[ind])
        pmax = adjpvalue[ind]
        # This null hypothesis has been processed
        processed[ind] = 1


        # Update the hypothesis weights after a null hypothesis has been processed
        temp = w
        for (j in 1:m) {
          if (processed[j] == 0)
            w[j] = temp[j] + temp[ind] * g[ind, j] else w[j] = 0
        }


        # Update the transition parameters (connection weights) after the rejection
        temp = g
        for (j in 1:m) {
          for (k in 1:m) {
            if (processed[j] == 0 & processed[k] == 0 & j != k & temp[j, ind] * temp[ind, j] != 1)
              g[j, k] = (temp[j, k] + temp[j, ind] * temp[ind, k])/(1 - temp[j, ind] * temp[ind, j]) else g[j, k] = 0
          }
        }
      }
      else {
        adjpvalue[which(processed==0)]=1
      }
    }
    result = adjpvalue
  }
  else if (call == TRUE) {
    w = paste0("Weight={",paste(round(par[[2]]$weight, 3), collapse = ","),"}")
    g = paste0("Transition matrix={",paste(as.vector(t(par[[2]]$transition)), collapse = ","),"}")
    result=list(list("Chain procedure"),list(w,g))
  }



  return(result)
}
# End of ChainAdj