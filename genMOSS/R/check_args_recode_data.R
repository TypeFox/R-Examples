check_args_recode_data <-
function (data, dimens, alpha) {

  if (!is.numeric(alpha)) {
    stop("class(alpha) != 'numeric'")
  }
  else if (length(alpha) != 1) {
    stop("length(alpha) != 1")
  }
  else if (alpha <= 0) {
    stop ("alpha must be positive")
  }
  if (!is.numeric(dimens)) {
    stop ("class(dimens) != 'numeric'")
  }
  else if (dim(data)[2] < 8) {
    stop ("dim(data)[2] must be at least 8")
  }
  else if (length(dimens) != dim(data)[2]) {
    stop ("dim(data)[2] != length(dimens)")
  }
  else if (dim(data)[1] < 2) {
    stop ("Excluding incomplete cases, dim(data)[1] must be at least 2")
  }
  if (dimens[length(dimens)] != 2) {
    stop ("Response must be binary")
  }
  for (i in 1:length(dimens)) {
    if (dimens[i] - floor(dimens[i]) != 0 || dimens[i] < 2) {
      stop("All entries of dimens must be integers greater than or equal to 2")
    }
  }

}
