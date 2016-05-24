######################################################################################################################

# Function: HochbergAdj.
# Argument: p, Vector of p-values (1 x m)
#           par, List of procedure parameters: vector of hypothesis weights (1 x m)
# Description: Hochberg multiple testing procedure.

HochbergAdj = function(p, par) {

  # Determine the function call, either to generate the p-value or to return description
  call = (par[[1]] == "Description")

  # Number of p-values
  m = length(p)

  # Extract the vector of hypothesis weights (1 x m)
  if (!any(is.na(par[[2]]))) {
    if (is.null(par[[2]]$weight)) stop("Analysis model: Hochberg procedure: Hypothesis weights must be specified.")
    w = par[[2]]$weight
  } else {
    w = rep(1/m, m)
  }

  if (any(call == FALSE) | any(is.na(call))) {
    # Error checks
    if (length(w) != m) stop("Analysis model: Hochberg procedure: Length of the weight vector must be equal to the number of hypotheses.")
    if (sum(w)!=1) stop("Analysis model: Hochberg procedure: Hypothesis weights must add up to 1.")
    if (any(w < 0)) stop("Analysis model: Hochberg procedure: Hypothesis weights must be greater than 0.")

    if (max(w) == min(w)){
      # Index of ordered pvalue
      ind <- order(p, decreasing = TRUE)

      # Adjusted p-values
      result <- pmin(1, cummin(cumsum(w[ind]) * p[ind]/w[ind]))[order(ind)]
    } else {

      # Compute the weighted incomplete Simes p-value for an intersection hypothesis
      incsimes<-function(p,u) {
        k<-length(u[u!=0 & !is.nan(u)])
        if (k>1) {
          temp=matrix(0,2,k)
          temp[1,]<-p[u!=0]
          temp[2,]<-u[u!=0]
          sort<-temp[,order(temp[1,])]
          modu<-u[u!=0]
          modu[1]<-0
          modu[2:k]<-sort[2,1:k-1]
          incsimes<-min((1-cumsum(modu))*sort[1,]/sort[2,])
        } else if (k==1)  {
          incsimes<-p[u!=0]/u[u!=0]
        } else if (k==0) incsimes<-1
        return(incsimes)
      }
      # End of incsimes

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
        int.pval[i, ] <- int[i, ] * incsimes(p, w.loc)
      }

      result <- apply(int.pval, 2, max)
    }
  }
  else if (call == TRUE) {
    weight = paste0("Weight={",paste(round(w,2), collapse = ","),"}")
    result=list(list("Hochberg procedure"),list(weight))
  }

  return(result)
}
# End of HochbergAdj