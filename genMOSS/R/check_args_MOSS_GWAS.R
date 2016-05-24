check_args_MOSS_GWAS <-
function (alpha, c, cPrime, q, replicates, maxVars, data, dimens, confVars, k) {

  if (!is.numeric(alpha)) {
    stop("class(alpha) != 'numeric'")
  }
  else if (length(alpha) != 1) {
    stop("length(alpha) != 1")
  }
  else if (alpha <= 0) {
    stop ("alpha must be positive")
  }
  else if (!is.numeric(c)) {
    stop ("class(c) != 'numeric'")
  }
  else if (length(c) != 1) {
    stop("length(c) != 1")
  }
  else if (!is.numeric(cPrime)) {
    stop ("class(cPrime) != 'numeric'")
  }
  else if (length(cPrime) != 1) {
    stop("length(cPrime) != 1")
  }
  else if (!is.numeric(q)) {
    stop ("class(q) != 'numeric'")
  }
  else if (length(q) != 1) {
    stop("length(q) != 1")
  }
  else if (!is.numeric(replicates)) {
    stop ("class(replicates) != 'numeric'")
  }
  else if (length(replicates) != 1) {
    stop("length(replicates) != 1")
  }
  else if (replicates - floor(replicates) != 0 || replicates <= 0) {
    stop ("Replicates must be a positive integer")
  }
  else if (!is.numeric(maxVars)) {
    stop ("class(maxVars) != 'numeric'")
  }
  else if (length(maxVars) != 1) {
    stop("length(maxVars) != 1")
  }
  else if (maxVars - floor(maxVars) != 0 || maxVars < 3 || maxVars > 6) {
    stop ("maxVars must be an integer between 3 and 6")
  }
  else if (!is.null(k) && !is.numeric(k)) {
    stop ("class(k) must be 'NULL' or 'numeric'")
  }
  else if (is.numeric(k) && length(k) != 1) {
    stop("k must be a single integer")
  }
  else if (is.numeric(k) && (k - floor(k) != 0 || k < 2)) {
    stop ("k must be 'NULL' or an integer greater than or equal to 2")
  }
  if (!is.numeric(dimens)) {
    stop ("class(dimens) != 'numeric'")
  }
  else if (dim(data)[2] < 8 || length(dimens) != dim(data)[2]) {
    stop ("dim(data)[2] must be at least 8 and be equal to length(dimens)")
  }
  else if (dim(data)[1] < 2) {
    stop ("Excluding incomplete cases, dim(data)[1] must be at least 2")
  }
  if (dimens[length(dimens)] != 2) {
    stop ("Response must be binary")
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
  else if(!is.null(confVars) && !is.character(confVars)) {
    stop ("class(confVars) must be 'NULL' or 'character'")
  }
  else if (length(confVars) < 0 || maxVars <= length(confVars) + 1) {
    stop ("length(confVars) must be between 0 and maxVars - 2")
  }
  else if (length(which(!(confVars %in% colnames(data)))) != 0) {
    stop(paste("Some confounding variables not found in colnames(data)"))
  }
  else if (colnames(data)[dim(data)[2]] %in% confVars) {
    stop(paste("The response cannot be in confVars")) 
  }
  for (i in 1:length(dimens)) {
    if (dimens[i] - floor(dimens[i]) != 0 || dimens[i] < 2) {
      stop("All entries of dimens must be integers greater than or equal to 2")
    }
  }
}
