formatp <- function(p, decimals = c(2, 3), cuts = 0.01, lowerbound = 0.001, 
                    leading0 = TRUE, avoid1 = FALSE) {
  
  # If any inputs are not correct class or out of range, return error
  if (min(p) < 0 | max(p) > 1) {
    stop("Please ensure that all p-values are between 0 and 1")
  }
  if (!is.numeric(decimals)) {
    stop("For decimals input, please enter numeric value or vector")
  }
  if (!is.numeric(cuts)) {
    stop("For cuts input, please enter numeric value or vector")
  }
  if (!is.numeric(lowerbound)) {
    stop("For lowerbound input, please enter numeric value")
  }
  if (!is.logical(leading0)) {
    stop("For leading0 input, please enter TRUE or FALSE")
  }
  if (!is.logical(avoid1)) {
    stop("For avoid1 input, please enter TRUE or FALSE")
  }
  
  # Determine number of decimal points to use
  if (length(decimals) > 1) {
    dec <- cut(x = p, breaks = c(0, cuts, 1), labels = rev(decimals), include.lowest = TRUE, right = FALSE)
  } else {
    dec <- rep(decimals, length(p))
  }
  
  # Create spf to control format of string
  spf <- paste("%0.", dec, "f", sep = "")
  
  # Format pval
  pval <- c()
  for (ii in 1:length(p)) {
    if (leading0 == TRUE) {
      if (p[ii] < lowerbound) {
        pval[ii] <- paste("<", as.numeric(lowerbound), sep = "")
      } else {
        pval[ii] <- sprintf(spf[ii], p[ii])
      }
      if (avoid1 == TRUE & unlist(strsplit(pval[ii], ""))[1] == "1") {
        pval[ii] <- paste(">0.", paste(rep("9", dec[ii]), sep = "", collapse = ""), sep = "")
      }
    } else {
      if (p[ii] < lowerbound) {
        pval[ii] <- paste("<", paste(strsplit(as.character(lowerbound), "")[[1]][-1], collapse = ""), sep = "")
      } else {
        first <- paste(strsplit(sprintf(spf[ii], p[ii]), "")[[1]])[1]
        if (first == "1") {
          pval[ii] <- paste(strsplit(sprintf(spf[ii], p[ii]), "")[[1]], collapse = "")
        } else {
          pval[ii] <- paste(strsplit(sprintf(spf[ii], p[ii]), "")[[1]][-1], collapse = "")
        }
      }
      if (avoid1 == TRUE & unlist(strsplit(pval[ii], ""))[1] == "1") {
        pval[ii] <- paste(">.", paste(rep("9", dec[ii]), sep = "", collapse = ""), sep = "")
      }
    }
  }
  
  # Return pval
  return(pval)
  
}