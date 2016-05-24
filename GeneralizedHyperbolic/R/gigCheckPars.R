### Check for boundary and error cases
### Return case normal, gamma, inverse gamma or error
### In case of error, include appropriate message
gigCheckPars <- function(param, ...) {

  param <- as.numeric(param)
  chi <- param[1]
  psi <- param[2]
  lambda <- param[3]

  case <- "normal"
  errMessage <- ""

  if (length(param) != 3) {
    case <- "error"
    errMessage <- "param vector must contain 3 values"
  } else {
    if (chi < 0) {
      case <- "error"
      errMessage <- "chi must be non-negative"
    } else if (psi < 0) {
      case <- "error"
      errMessage <- "psi must be non-negative"
    } else if (isTRUE(all.equal(0, chi, ...))) {
      if (lambda <= 0) {
        case <- "error"
        errMessage <- "lambda must be positive when chi = 0"
      } else if (isTRUE(all.equal(0, psi, ...))) {
        case <- "error"
        errMessage <- "psi and chi cannot both be 0"
      } else {
        case <- "gamma"
      }
    } else if (isTRUE(all.equal(0, psi, ...))) {
      if (lambda >= 0) {
        case <- "error"
        errMessage <- "lambda must be negative when psi = 0"
      } else if (isTRUE(all.equal(0, chi, ...))) {
        case <- "error"
        errMessage <- "psi and chi cannot both be 0"
      } else {
        case <- "invgamma"
      }
    }
  }

  result <- list(case = case, errMessage = errMessage)
  return(result)
}
