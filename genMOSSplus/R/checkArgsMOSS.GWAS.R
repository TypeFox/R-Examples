checkArgsMOSS.GWAS <-
function (alpha, c, cPrime, q, replicates, data, maxVars, dimens, k, seed) {

  if (!is.numeric(alpha)) {
    stop("class(alpha) != 'numeric'")
  }
  else if (alpha <= 0) {
    stop ("alpha must be positive")
  }
  else if (!is.numeric(c)) {
    stop ("class(c) != 'numeric'")
  }
  else if (!is.numeric(cPrime)) {
    stop ("class(c) != 'numeric'")
  }
  else if (!is.numeric(q)) {
    stop ("class(c) != 'numeric'")
  }
  else if (!is.numeric(replicates)) {
    stop ("class(replicates) != 'numeric'")
  }
  else if (replicates - floor(replicates) != 0 || replicates <= 0) {
    stop ("replicates must be a positive integer")
  }
  else if (!is.numeric(maxVars)) {
    stop ("class(maxVars) != 'numeric'")
  }
  else if (maxVars - floor(maxVars) != 0 || maxVars < 3 || maxVars > 6) {
    stop ("maxVars must be an integer between 3 and 6")
  }
  else if (!is.null(k) && !is.numeric(k)) {
    stop ("class(k) must be 'NULL' or 'numeric'")
  }
  else if (is.numeric(k) && (k - floor(k) != 0 || k < 2)) {
    stop ("k must be 'NULL' or an integer greater than or equal to 2")
  }
  else if (dim(data)[1] < 2) {
    stop ("Excluding incomplete cases, dim(data)[1] must be at least 2")
  }
  else if (dim(data)[2] < 3) {
    stop ("dim(data)[2] must be at least 3")
  }
  else if (c < 0 || c > 1) {
    stop ("c must be between 0 and 1")
  }
  else if (cPrime > c || cPrime < 0) {
    stop ("cPrime must be between 0 and c")
  }
  else if (q < 0 || q > 1) {
    stop ("q must be between 0 and 1")
  }  
  if (!is.numeric(dimens)) {
    stop ("class(dimens) != 'numeric'")
  }
  if (dimens[length(dimens)] != 2) {
    stop ("response must be binary")
  }
  for (i in 1:length(dimens)) {
    if (dimens[i] - floor(dimens[i]) != 0 || dimens[i] < 2) {
      stop("all entries of dimens must be integers greater than or equal to 2")
    }
  }
}
