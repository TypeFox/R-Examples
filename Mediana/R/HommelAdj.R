######################################################################################################################

# Function: HommelAdj.
# Argument: p, Vector of p-values (1 x m)
#           par, List of procedure parameters: vector of hypothesis weights (1 x m)
# Description: Hommel multiple testing procedure (using the closed testing procedure).

HommelAdj = function(p, par) {

  # Determine the function call, either to generate the p-value or to return description
  call = (par[[1]] == "Description")


  # Number of p-values
  m = length(p)

  # Extract the vector of hypothesis weights (1 x m)
  if (!any(is.na(par[[2]]))) {
    if (is.null(par[[2]]$weight)) stop("Analysis model: Hommel procedure: Hypothesis weights must be specified.")
    w = par[[2]]$weight
  } else {
    w = rep(1/m, m)
  }

  if (any(call == FALSE) | any(is.na(call))) {
    # Error checks
    if (length(w) != m) stop("Analysis model: Hommel procedure: Length of the weight vector must be equal to the number of hypotheses.")
    if (sum(w)!=1) stop("Analysis model: Hommel procedure: Hypothesis weights must add up to 1.")
    if (any(w < 0)) stop("Analysis model: Hommel procedure: Hypothesis weights must be greater than 0.")

    # Weighted Simes p-value for intersection hypothesis
    simes <- function(p, w) {
      nb <- length(w[w != 0 & !is.nan(w)])
      if (nb > 1) {
        p.sort <- sort(p[w != 0])
        w.sort <- w[w != 0][order(p[w != 0])]
        simes <- min(p.sort/cumsum(w.sort))
      }
      else if (nb == 1)
        simes <- pmin(1, p/w)
      else if (nb == 0)
        simes <- 1
      return(simes)
    }

    # number of intersection
    nbint <- 2^m - 1

    # matrix of intersection hypotheses
    int <- matrix(0, nbint, m)
    for (i in 1:m) {
      for (j in 0:(nbint - 1)) {
        k <- floor(j/2^(m - i))
        if (k/2 == floor(k/2))
          int[j + 1, i] <- 1
      }
    }

    # matrix of local p-values
    int.pval <- matrix(0, nbint, m)

    # vector of weights for local test
    w.loc <- rep(0, m)

    # local p-values for intersection hypotheses
    for (i in 1:nbint) {
      w.loc <- w * int[i, ]/sum(w * int[i, ])
      int.pval[i, ] <- int[i, ] * simes(p, w.loc)
    }

    adjpvalue <- apply(int.pval, 2, max)
    result = adjpvalue
  }
  else if (call == TRUE) {
    weight = paste0("Weight={",paste(round(w,2), collapse = ","),"}")
    result=list(list("Hommel procedure"),list(weight))
  }

  return(result)
}
# End of HommelAdj